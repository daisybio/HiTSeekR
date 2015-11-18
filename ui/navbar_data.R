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
  do.call(wellPanel, list(fluidRow(column(6, 
                                          selectInput("dataSelectedNormalization", "Normalization", normalizationChoices(), "Raw")
  ), 
  column(6, 
         checkboxInput("dataShowHelpText", "Show help text", FALSE)
  )
  )
  ))
})

output$uiOutput_data <- renderUI({ 
  exp.data <- processedData()
  plates <- unique(as.character(exp.data$Plate))
  replicates <- unique(as.character(exp.data$Replicate))
  elements <- list(  
    tabPanel("Normalized Data", 
             conditionalPanel("input.dataShowHelpText",
                              HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                         normDataTableInfoText,
                                         '</div>', sep="")
                              )
             ), dataTableOutput("table_processedData")),
    tabPanel("Plate Signal Variation",
             conditionalPanel("input.dataShowHelpText",
                              HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                         plateMeanInfoText,
                                         '</div>', sep="")
                              )
             ), plotOutput("plateMeanPlot2", height="auto")),
    tabPanel("Whole Screen Scatter Plot", 
             conditionalPanel("input.dataShowHelpText",
                              HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                         wholeScreenInfoText,
                                         '</div>', sep="")
                              )
             ), uiOutput("scatterPlotTagList")),
    tabPanel("Signal Distribution", 
             conditionalPanel("input.dataShowHelpText",
                              HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                         signalDistributionInfoText,
                                         '</div>', sep="")
                              )
             ), plotOutput("signalDistPlot", height="auto")),
    tabPanel("QQ Plot", 
             conditionalPanel("input.dataShowHelpText",
                              HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                         qqPlotInfoText,
                                         '</div>', sep="")
                              )
             ), plotOutput("signalqqPlot", height="auto")),
    tabPanel("Plate Viewer",
             conditionalPanel("input.dataShowHelpText",
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
               tags$script(src = "https://code.highcharts.com/modules/heatmap.js"),
               showOutput("intHeatmapPlot", "highcharts"),
               showOutput("intPlateScatterPlot", "highcharts")
             )
    )
    
  )
  
  if(!is.null(input$replicateCol)) 
  {
    elements <- c(elements, list(
      tabPanel("Replicate Correlation",              
               conditionalPanel("input.dataShowHelpText",
                                HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                           replicateCorrInfoText,
                                           '</div>', sep="")
                                )
               ), plotOutput("replicateCorrPlot2"))
    ))
  }
  do.call(tabsetPanel, elements)
})