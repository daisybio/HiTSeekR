output$scatterPlotTagList <- renderUI({
  exp.data <- processedData()
  if(is.null(exp.data)) return(NULL)
  if(nrow(exp.data) > 5000){
    plots <- plotOutput("scatterPlot")
  } 
  else{
    plots <- foreach(experiment = unique(exp.data$Experiment)) %do% {
      plot.name <- paste(experiment, "IntScatterPlot", sep="")
      showOutput(plot.name, "highcharts")
    }
    do.call(tagList, plots)
  }
})

output$uiOutput_data <- renderUI({ 
elements <- list(
  tabPanel("Raw Data", dataTableOutput("table_rawData")),
  tabPanel("Processed Data", dataTableOutput("table_processedData")),
  tabPanel("Scatter Plots",  uiOutput("scatterPlotTagList"))
)

do.call(tabsetPanel, elements)
})