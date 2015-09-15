output$uiOutput_quality_control <- renderUI({  
  
  elements <- list(
    tabPanel("Plate Signal Spread",
      plotOutput("plateMeanPlot"),
      shinyalert("plateMeanInfo")
    ),
    tabPanel("Row and Column Effect",
      plotOutput("rowAndColumn"),
      shinyalert("rowColumnInfo")
    )
  )
  
  if(input$hasControls){
    elements <- c(elements, list(
    tabPanel("Control Signal Spread", 
      plotOutput("controlPlot")
    ),
    tabPanel("Control Separability",
      plotOutput("controlPerformancePlot")
    )))
  }
  if(!is.null(input$replicateCol)) 
  {
    elements <- c(elements, list(
      tabPanel("Replicate Correlation",
               plotOutput("replicateCorrPlot")
    )))
  }
  
  do.call(tabsetPanel, elements)
})
