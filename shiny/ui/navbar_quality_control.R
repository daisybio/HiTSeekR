output$uiOutput_quality_control <- renderUI({  
  
  elements <- list(
    tabPanel("Plate Signal Variation",
      plotOutput("plateMeanPlot"),
      conditionalPanel("input.showHelpPages",
                       HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                  plateMeanInfoText,
                                  '</div>', sep="")
                       )
      )
    ),
    tabPanel("Row and Column Effect",
      plotOutput("rowAndColumn"),
      conditionalPanel("input.showHelpPages",
                       HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                  rowColumnInfoText,
                                  '</div>', sep="")
                       )
      )    )
  )
  
  if(input$hasControls){
    elements <- c(elements, list(
    tabPanel("Control Signal Spread", 
      plotOutput("controlPlot"),
      conditionalPanel("input.showHelpPages",
                       HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                  controlPlotInfoText,
                                  '</div>', sep="")
                       )
      )
    ),
    tabPanel("Control Separability",
      plotOutput("controlPerformancePlot"),
      conditionalPanel("input.showHelpPages",
                       HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                  controlPerformancePlotInfoText,
                                  '</div>', sep="")
                       )
      )
    )))
  }
  if(!is.null(input$replicateCol)) 
  {
    elements <- c(elements, list(
      tabPanel("Replicate Correlation",
               plotOutput("replicateCorrPlot"),
               conditionalPanel("input.showHelpPages",
                                HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                           replicateCorrInfoText,
                                           '</div>', sep="")
                                )
               )
    )))
  }
  
  do.call(tabsetPanel, elements)
})
