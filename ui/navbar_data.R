output$scatterPlotTagList <- renderUI({
  exp.data <- processedData()
  if(is.null(exp.data)) return(NULL)
  if(nrow(exp.data) > 5000){
    plots <- plotOutput("scatterPlot")
  } 
  else{
    plots <- foreach(experiment = unique(exp.data$Experiment)) %do% {
      foreach(readout = unique(exp.data$Readout)) %do% {
        plot.name <- paste(experiment, readout, "IntScatterPlot", sep="")
        showOutput(plot.name, "highcharts")
      }
    }
    do.call(tagList, plots)
  }
})

output$uiOutput_dataWellPanel <- renderUI({
  slt <- selectInput("dataSelectedNormalization", "Normalization", normalizationChoices(), "Raw")
  do.call(wellPanel, list(slt))
})

output$uiOutput_data <- renderUI({ 
  exp.data <- processedData()
  plates <- unique(as.character(exp.data$Plate))
  replicates <- unique(as.character(exp.data$Replicate))
elements <- list(  
  tabPanel("Normalized Data", dataTableOutput("table_processedData")),
  tabPanel("Plate Signal Variation", plotOutput("plateMeanPlot2")),
  tabPanel("Replicate Correlation", plotOutput("replicateCorrPlot2")),
  tabPanel("Whole Screen Scatter Plot",  uiOutput("scatterPlotTagList")),
  tabPanel("Signal Distribution", plotOutput("signalDistPlot")),
  tabPanel("QQ Plot", plotOutput("signalqqPlot")),
  tabPanel("Plate Viewer",
           sidebarPanel(
             selectInput("plateSelected", "Select a plate:", plates, plates[1]),
             selectInput("replicateSelected", "Select a replicate (heatmap):", replicates, replicates[1]),
             selectInput("heatmapExperimentSelected", "Select an experiment:", experiments(), experiments()[1]),
             selectInput("heatmapReadoutSelected", "Select a readout:", readouts(), readouts()[1])             
           ),mainPanel(
             tags$script(src = "https://code.highcharts.com/modules/heatmap.js"),
             showOutput("intHeatmapPlot", "highcharts"),
             showOutput("intPlateScatterPlot", "highcharts")
           )
  )
  
)

do.call(tabsetPanel, elements)
})