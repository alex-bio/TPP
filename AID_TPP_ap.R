##Analysis of Independent Differences (AID)
##For TPP data. Written by Alexandra Panov.
##Written in two steps: first is AID, a function of control and treatment files, peptide_threshold, pvalue_threshold, and a vector of temperature values. Second is the Shiny App.
##Required packages.
list.of.packages <- c("shiny", "writexl","ggplot2","Cairo","MASS","mvtnorm","plyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(shiny)
library(writexl)
library(ggplot2)
library(Cairo)
library(MASS)
library(mvtnorm)
library(plyr)

AID <- function(control, treatment, peptide_threshold,pvalue_threshold,temp_vec) {
  ####Filtering steps.#### 
  ##First, get rid of reverse hits and NAs. Any row with at least one NA value will be eliminated.
  control_filtered <- control[!grepl("##", control$Protein.Id),]
  control_filtered <- control_filtered[complete.cases(control_filtered),]
  treatment_filtered <- treatment[!grepl("##", treatment$Protein.Id),]
  treatment_filtered <- treatment_filtered[complete.cases(treatment_filtered),]
  
  ##Get rid of all proteins that have peptide values < 3. 3 is default, but can change if you want. Note that this number is inclusive. i.e. Peptides > or = are passed forward and only < are eliminated. This step will be eliminated when plugging into the Shiny App because it will be user defined.
  #peptide_threshold <- 3
  
  ##Get rid of any observations with less peptides than set by the peptide threshold.
  control_filtered <- subset(control_filtered, control_filtered$Number.of.peptides >= peptide_threshold)
  treatment_filtered <- subset(treatment_filtered, treatment_filtered$Number.of.peptides >= peptide_threshold)
  
  ##Get rid of duplicates. If the protein ID is, for example, a gene symbol, multiple isoforms may be present with the same gene symbol name/string. This analysis cannot tell the difference, so any non-unique observations are removed.
  control_filtered <- control_filtered[!(duplicated(control_filtered$Protein.Id) | duplicated(control_filtered$Protein.Id, fromLast = TRUE)), ]
  treatment_filtered <- treatment_filtered[!(duplicated(treatment_filtered$Protein.Id) | duplicated(treatment_filtered$Protein.Id, fromLast = TRUE)), ]
  
  ##Merge the data by protein ID, then split again later to have individual control and treatment dataframes. Reordering ensures that control and treatment dataframes match Protein.Id by row.
  all_data <- merge(control_filtered,treatment_filtered, by="Protein.Id")
  in_common <- data.frame(Protein.Id = all_data[,c("Protein.Id")])
  control_filtered <- merge(in_common,control_filtered, by="Protein.Id")
  control_filtered <- control_filtered[order(control_filtered$Protein.Id),]
  treatment_filtered <- merge(in_common,treatment_filtered, by="Protein.Id")
  treatment_filtered <- treatment_filtered[order(treatment_filtered$Protein.Id),]
  
  ##Pick out the columns with values for fraction non-denatured protein. All columns with "temp" are selected.
  control_temps <- control_filtered[, grepl("temp", names(control_filtered))]
  treatment_temps <- treatment_filtered[, grepl("temp", names(treatment_filtered))]
  
  ##Find the max in each row and normalize, setting the max value = 1.
  control_temps_max <- apply(control_temps,1,max) 
  control_norm <- control_temps/control_temps_max
  treatment_temps_max <- apply(treatment_temps,1,max)
  treatment_norm <- treatment_temps/treatment_temps_max
  
  ##If the user column order doesn't match between control and treatment, this step will fix it.
  control_norm <- control_norm[,order(colnames(control_norm),decreasing=FALSE)]
  treatment_norm <- treatment_norm[,order(colnames(treatment_norm), decreasing=FALSE)]
  ##Now subtract!
  deltas <- treatment_norm - control_norm
  
  ##Now pick the middle 95% of deltas (differences between conditions of fraction non-denatured protein) for each temperature channel and calculate the Normal fits. SDs will come up later when calculating individual Normal p-value (pnorm) for each temperature channel.
  middle_data <- vector("list",ncol(deltas))
  means_from_fit <- c(1:ncol(deltas))
  sds_from_fit <- c(1:ncol(deltas))
  
  for (i in 1:ncol(deltas)) {
    quant_low <- quantile(deltas[,i], 0.025)
    quant_high <- quantile(deltas[,i], 0.975)
    middle_data[[i]] <- (subset(deltas[,i], deltas[,i] > quant_low & deltas[,i] < quant_high))
    x <- middle_data[[i]]
    fit <- fitdistr(x,"normal")
    means_from_fit[i] <- fit[[1]][[1]]
    sds_from_fit[i] <- fit[[1]][[2]]
  }
  ##Only need the middle data for calculating covariance. Use this covariance because we used the middle 95% observations to calculate the Normal parameters for each temperature channel.
  middle_deltas <- as.data.frame(do.call(cbind,middle_data))
  cov_middle <- cov(middle_deltas)
  
  ##Set up the integration bounds. Want area from tail inwards, so have to ask if each delta/difference value is less than the mean for that temperature channel.
  lower_bounds <- vector("list",nrow(deltas))
  upper_bounds <- vector("list",nrow(deltas))
  for (i in 1:nrow(deltas)) {
    lower_bounds[[i]] <- ifelse(deltas[i,] < means_from_fit, -Inf, deltas[i,])
    upper_bounds[[i]] <- unlist(ifelse(deltas[i,] < means_from_fit, deltas[i,], Inf))
  }
  ##Applying pmvnorm, a function written to calculate the Multivariate Normal p-value using four arguments: a vector of lower bounds, a vector of upper bounds, a vector of means, and the covariance matrix (can use correlation or covariance). The covariance is from the middle 95% of the data, in accordance with the Multivariate Normal fit. Each argument vector is of length = number of temperature channels.
  pvals <- vector("list",nrow(deltas))
  lower <- matrix(unlist(lower_bounds),ncol=ncol(deltas),byrow=T)
  upper <- matrix(unlist(upper_bounds),ncol=ncol(deltas),byrow=T)
  
  for (i in 1:nrow(deltas)) {
    lower_b <- lower[i,]
    upper_b <- upper[i,]
    pvals[[i]] <- pmvnorm(lower=lower_b, upper=upper_b, mean=means_from_fit, sigma=cov_middle)
  }
  pvals_numeric <- unlist(pvals)
  ##Multiply by 2^(ncol(deltas)) because each value in the Multivariate Observation is equally likely to be stabilized or destabilized. Can think of as: equally likely to be as far into the tail on one side of the distribution versus the other side. This is explained again below for individual p-value calcuations. For Multivariate Normal p-values that equal -Inf, the value is listed as -1000. Note that this step does not alter the ranking/p-value order.
  pvals_numeric <- pvals_numeric*(2^(ncol(deltas)))
  log_pvals <- log(pvals_numeric)
  log_pvals <- ifelse(log_pvals ==-Inf, -1000,log_pvals)
  
  indiv_pvals <- vector("list", ncol(deltas))
  ##Calculating individual p-values. Note that we are interested in how far the delta/difference is into the tail. If the pnorm p-value calculated is > 0.5, then subtract the p-value from 1 because we want the integral from the value of the delta to +Inf. We multiply by 2 because the delta value is equally likely to be either stabilized or destabilized.
  for (i in 1:ncol(deltas)){
    x <- deltas[,i]
    y <- pnorm(x,mean=means_from_fit[i],sd=sds_from_fit[i])
    indiv_pvals[[i]] <- 2*ifelse(y < 0.5,y,1-y)
  }
  all_indiv_pvals <- data.frame(t(matrix(unlist(indiv_pvals),ncol=nrow(deltas),byrow=T)))
  colnames(all_indiv_pvals) <- paste("pval_temp", formatC(1:ncol(all_indiv_pvals), width=nchar(ncol(deltas)), flag="0"), sep="_")
  ##Counting the number of temperature channels for each protein that have an individual p-value <= pvalue_threshold. Defining the pvalue_threshold is commented out because it is user defined in the Shiny app.
  #pvalue_threshold <- 0.05
  pval_count <- rowSums(all_indiv_pvals <= pvalue_threshold)
  ##Issuing a variance warning if variance of each protein observation is more than 2 standard deviations away from the mean.
  var_deltas <- apply(deltas,1,var)
  var_warning <- ifelse(var_deltas > (mean(var_deltas)+(2*sd(var_deltas))), "Warning: high variance", ifelse(var_deltas < (mean(var_deltas)-(2*sd(var_deltas))), "Warning: low variance", ""))
  ##Summing the signs
  sum_of_signs <- apply(deltas, 1, function(x) sum(sign(x)))
  
  ##Now we can build the results.
  deltas_final <- cbind(Protein.Id = control_filtered$Protein.Id,deltas,log_Multiv_Norm_pval=log_pvals,all_indiv_pvals,pval_count=pval_count,var_warning=var_warning,sum_of_signs=sum_of_signs)
  deltas_final$prediction <- ifelse(deltas_final$sum_of_signs < 0, "destabilized", ifelse(deltas_final$sum_of_signs > 0, "stabilized", "undetermined"))
  ##Sort first by log of Multivariate Normal p-value, from least to most, then by pval_count, or number of individual Normal p-values that are below the user defined p-value threshold, and then by magnitude of sum_of_signs, which indicates approximate directionality.
  deltas_final <- deltas_final[order(deltas_final$log_Multiv_Norm_pval,-deltas_final$pval_count,-abs(deltas_final$sum_of_signs)),]
  
  ##Making the Normalized Data spreadsheet.
  unnorm_data <- all_data
  unnorm_data <- unnorm_data[order(unnorm_data$Protein.Id),]
  colnames(unnorm_data) <- gsub(".x","_control",colnames(unnorm_data))
  colnames(unnorm_data) <- gsub(".y","_treatment",colnames(unnorm_data))
  control_temp <- control_norm
  colnames(control_temp) <- paste(colnames(control_temp), "control", sep = "_")
  treatment_temp <- treatment_norm
  colnames(treatment_temp) <- paste(colnames(treatment_temp), "treatment", sep = "_")
  norm_data2 <- cbind(control_temp,treatment_temp)
  half <- unnorm_data[ ,!(names(unnorm_data) %in% colnames(norm_data2))]
  normalized_data <- cbind(half,norm_data2)
  
  ##This scratchwork below makes a dataframe that shiny can use to graph the data.
  scratchwork <- t(norm_data2)
  colnames(scratchwork) <- normalized_data$Protein.Id
  ##The Temperature vector, temp_vec, is user defined in the Shiny app. These are the actual temperatures used in the experiment. These values are only used for graphing purposes.
  #temp_vec <- c(37,40,42,47,50,52,56,58,60,64)
  ready_for_shiny <- data.frame(scratchwork)
  ready_for_shiny$Temperature <- rep(temp_vec,2)
  N <- (nrow(ready_for_shiny)/2)
  ready_for_shiny <- cbind(tpp.cond=c(rep("control",N),rep("treatment",N)), ready_for_shiny)

return(list("AID_output"=deltas_final,"normalized_data"=normalized_data,"ready_for_shiny"=ready_for_shiny))
}
##The write_xlsx step is commented out because this also happens inside the Shiny app.
#write_xlsx(list(AID_output = deltas_final, normalized_data = normalized_data), "./I_heart_MS_4ever.xlsx")
##End of AID.

##Start of Shiny App.
ui <- shinyUI(fluidPage(
  titlePanel("AID: a TPP Data Analysis"),
  sidebarLayout(
    sidebarPanel(
      helpText("Analysis of Independent Differences (AID) for analyzing thermal proteome profiling data. For more information, please see Panov, A. & Gygi, S.P. (2019)."),
      fileInput("control", "Choose Control CSV File", accept=c("text/csv","text/comma-separated-values,text/plain",".csv")),
      fileInput("treatment", "Choose Treatment CSV File", accept=c("text/csv","text/comma-separated-values,text/plain",".csv"))
      ,
      numericInput("peptide_threshold", label="Peptide Threshold",value=3, min=1,step=1),
      helpText("Peptide threshold is inclusive. If you pick 3, then all proteins with 3 or greater peptides are included."),
      numericInput("pvalue_threshold",label="P-value Threshold",value=0.05,min=0,max=1),
      helpText("While the ranking of shifted proteins incorporates this value, all proteins are shown in the analysis, regardless of this value."),
      textInput('temp_vec',"Enter temperatures separated by commas", "37,40,42,47,50,52,56,58,60,64"),
      hr(),
      conditionalPanel("output.files_ready",
                       actionButton("run_analysis",label="Run AID")),
      hr()),
    mainPanel(
      uiOutput("downloader"),
      uiOutput("protein_select"),
      uiOutput("plot_protein")
    )
   )
  )
)

server <- function(input,output,session){
  options(shiny.usecairo=T)
  output$files_ready <- reactive({
    return(!is.null(input$control) & !is.null(input$treatment))
  })
  outputOptions(output,"files_ready", suspendWhenHidden=FALSE)
  tpp.df <- NULL
  tpp.df <- eventReactive(input$run_analysis, {
    withProgress(message='Calculating',
                 detail='May take a few minutes...',value=0,{
                   if (is.null(input$control))
                     return(NULL)
                   if (is.null(input$treatment))
                     return(NULL)
                   control <- read.csv(input$control$datapath, header=TRUE)
                   incProgress(1/6)
                   treatment <- read.csv(input$treatment$datapath, header=TRUE)
                   incProgress(1/6)
                   temp_vec <- as.numeric(unlist(strsplit(input$temp_vec,",")))
                   AID(control,treatment,input$peptide_threshold,input$pvalue_threshold,temp_vec)
                 }
    )
  })
  
  output$downloader <- renderUI({
    if(!is.null(tpp.df())) 
      downloadButton('OutputFile', "Download Analysis File")
  })
  output$OutputFile <- downloadHandler(
    filename=function() {
      paste("TPP_analysis","_",Sys.time(),".xlsx",sep="")
    },
    content=function(file){
      write_xlsx(tpp.df()[c("AID_output","normalized_data")],file)
    }
  )
  output$protein_select <- renderUI({
    if(!is.null(tpp.df()))
  selectizeInput("protein","Protein:",
                 choices=colnames(within(tpp.df()[["ready_for_shiny"]],rm('tpp.cond','Temperature'))),
                 options=list(
                 placeholder="Type gene symbol",
                 maxOptions = 20000),
                 selected=NULL,
                 multiple=FALSE)
  #updateSelectizeInput(session,"protein",choices=colnames(within(tpp.df()[["ready_for_shiny"]],rm('tpp.cond','Temperature'))), server=TRUE,selected=NULL)
  })
  output$plot_protein <- renderUI({
    if(!is.null(tpp.df()))
      output$actual_plot <- renderPlot({
        p=  ggplot() +
          geom_line(data=tpp.df()[["ready_for_shiny"]], aes(x=Temperature,y=tpp.df()[["ready_for_shiny"]][,input$protein],group=tpp.cond,color=tpp.cond),size=1.5) +
          geom_point(data=tpp.df()[["ready_for_shiny"]], aes(x=Temperature,y=tpp.df()[["ready_for_shiny"]][,input$protein],group=tpp.cond,color=tpp.cond),size=4.5) +
          xlab('Temperature') +
          ylab('Fraction Non-denatured') +
          ylim(0,1) +
          ggtitle(input$protein)+
          scale_color_manual(values=c("#0072B2","#D55E00"), name="",labels=c("Control","Treatment")) +
          theme_bw()+
          theme(axis.title.x = element_text(size=17), axis.title.y = element_text(size=17), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.text = element_text(size=15),plot.title = element_text(size = 20, face = "bold"))
        print(p)
      })
    plotOutput("actual_plot")
  })
 
}
  
shinyApp(ui,server)

##End of Shiny App.