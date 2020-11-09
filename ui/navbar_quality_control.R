output$uiOutput_quality_control <- renderUI({  
  
  elements <- list(
    tabPanel("Plate Signal Variation",
       conditionalPanel("input.showHelpText",
                        HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                   plateMeanInfoText,
                                   '</div>', sep="")
                        )
       ),
      plotOutput("plateMeanPlot", height="auto")
    ),
    tabPanel("Row and Column Effect",
      conditionalPanel("input.showHelpText",
                       HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                  rowColumnInfoText,
                                  '</div>', sep="")
                       )
      ),
      plotOutput("rowAndColumn", height="auto")
    )
  )
  
  if(input$hasControls){
    elements <- c(elements, list(
    tabPanel("Control Signal Spread", 
      conditionalPanel("input.showHelpText",
                       HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                  controlPlotInfoText,
                                  '</div>', sep="")
                       )
      ),
      plotOutput("controlPlot", height="auto")
    ),
    tabPanel("Control Separability",
      checkboxInput("zfactor", label = "Use Z-factor instead of SSMD", FALSE),
      conditionalPanel("!input.zfactor",
                       conditionalPanel("input.showHelpText",
                                        HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                                   controlPerformancePlotInfoText,
                                                   '</div>', sep="")
                           )
                      ),
                      plotOutput("controlPerformancePlot", height="auto")
      ),
      conditionalPanel("input.zfactor",
                       conditionalPanel("input.showHelpText",
                                        HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                                   zScoreControlPerformancePlotInfoText,
                                                   '</div>', sep="")
                                        )
                       ),
                       plotOutput("zFactorControlPerformancePlot", height="auto")
      )
    )))
  }
  if(!is.null(input$replicateCol)) 
  {
    elements <- c(elements, list(
      tabPanel("Replicate Correlation",
               conditionalPanel("input.showHelpText",
                                HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                           replicateCorrInfoText,
                                           '</div>', sep="")
                                )
               ),
               checkboxInput("show.replicate.correlation", label = "Show correlation in plot", TRUE),
               plotOutput("replicateCorrPlot", height="auto")
    )))
  }
  
  do.call(tabsetPanel, elements)
})
