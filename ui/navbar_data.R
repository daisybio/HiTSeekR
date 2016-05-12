output$scatterPlotTagList <- renderUI({
  exp.data <- processedData()
  if(is.null(exp.data)) return(NULL)
  #if(nrow(exp.data) > 5000){
  plotOutput("scatterPlot")
  #   } 
  #   else{
  #     plots <- foreach(experiment = unique(exp.data$Experiment)) %do% {
  #       foreach(readout = unique(exp.data$Readout)) %do% {
  #         plot.name <- paste(experiment, readout, "IntScatterPlot", sep="")
  #         showOutput(plot.name, "highcharts")
  #       }
  #     }
  #     do.call(tagList, plots)
  #   }
})

output$uiOutput_dataWellPanel <- renderUI({
  do.call(wellPanel, list(selectInput("dataSelectedNormalization", "Normalization", normalizationChoices(), "Raw")))
})

output$uiOutput_data <- renderUI({ 
  exp.data <- processedData()
  plates <- unique(as.character(exp.data$Plate))
  replicates <- unique(as.character(exp.data$Replicate))
  elements <- list(  
    tabPanel("Normalized Data", 
             conditionalPanel("input.showHelpText",
                              HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                         normDataTableInfoText,
                                         '</div>', sep="")
                              )
             ), dataTableOutput("table_processedData")),
    tabPanel("Plate Signal Variation",
             conditionalPanel("input.showHelpText",
                              HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                         plateMeanInfoText,
                                         '</div>', sep="")
                              )
             ), plotOutput("plateMeanPlot2", height="auto")),
    tabPanel("Whole Screen Scatter Plot", 
             conditionalPanel("input.showHelpText",
                              HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                         wholeScreenInfoText,
                                         '</div>', sep="")
                              )
             ), uiOutput("scatterPlotTagList")),
    tabPanel("Signal Distribution", 
             conditionalPanel("input.showHelpText",
                              HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                         signalDistributionInfoText,
                                         '</div>', sep="")
                              )
             ), plotOutput("signalDistPlot", height="auto")),
    tabPanel("QQ Plot", 
             conditionalPanel("input.showHelpText",
                              HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                         qqPlotInfoText,
                                         '</div>', sep="")
                              )
             ), plotOutput("signalqqPlot", height="auto")),
    tabPanel("Plate Viewer",
             conditionalPanel("input.showHelpText",
                              HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                         plateViewerInfoText,
                                         '</div>', sep="")
                              )
             ), 
             sidebarPanel(
               selectInput("plateSelected", "Select a plate:", plates, plates[1]),
               selectInput("replicateSelected", "Select a replicate (heatmap):", replicates, replicates[1]),
               selectInput("heatmapExperimentSelected", "Select an experiment:", experiments(), experiments()[1]),
               selectInput("heatmapReadoutSelected", "Select a readout:", readouts(), readouts()[1])             
             ),mainPanel(
               tags$script(src = "http://code.highcharts.com/4.2.3/modules/heatmap.js"),
               showOutput("intHeatmapPlot", "highcharts"),
               showOutput("intPlateScatterPlot", "highcharts")
             )
    )
    
  )
  
  if(!is.null(input$replicateCol)) 
  {
    elements <- c(elements, list(
      tabPanel("Replicate Correlation",              
               conditionalPanel("input.showHelpText",
                                HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                           replicateCorrInfoText,
                                           '</div>', sep="")
                                )
               ), plotOutput("replicateCorrPlot2"))
    ))
  }
  do.call(tabsetPanel, elements)
})