KPM.run <- reactive({
  indicator.matrix <- targets.indicator.matrix()  
  result <- call.KPM(list(indicator.matrix))
  return(result)
})

#output$KPM.graph <- 

output$KPM.test <- renderPrint({
  input$startKPMButton
  
  result <- isolate(KPM.run())
  print(result)
})