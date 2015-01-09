#identifier column types
accTypes <- reactive({
  #if(input$screenType == "siRNA") return(c("RefSeq", "Entrez", "GeneSymbol"))
  #else if(input$screenType == "miRNA") return(c("MI", "MIMAT"))
  return(c("RefSeq", "Entrez", "MIMAT", "MI"))
})

ctrlTypes <- reactive({
  data <- rawData()
  if(!is.null(input$controlCol)) ctrlCol <- data[,input$controlCol]
  
  else ctrlCol <- data[,dataOptionDefaults()[["controlCol"]]]
  
  unique(as.character(ctrlCol))
})

negCtrl <- reactive({
  if(is.null(input$negCtrl)){
    return(dataOptionDefaults()[["negCtrls"]])
  }
})

posCtrl <- reactive({
  if(is.null(input$posCtrl)){
    return(dataOptionDefaults()[["posCtrls"]])
  }
})

output$uiOutput_data_options <- renderUI({
  elements <- list(column(4,
                          selectInput("screenType", "Type of screen", c("Gene knockout (e.g. siRNA)" = "siRNA", "miRNA inhibitor / mimics" = "miRNA"), dataOptionDefaults()[["screenType"]]),
                          selectInput("sampleCol", "Sample Name Column", dataColumns(), dataOptionDefaults()[["sampleCol"]]),
                          selectInput("plateCol", "Plate Column", dataColumns(), dataOptionDefaults()[["plateCol"]]),
                          selectInput("positionColType", "Position Column Type", c("Alpha well names" = "alpha", "Numeric" = "numeric"), dataOptionDefaults()[["posColType"]]),
                          selectInput("positionCol", "Position Column", dataColumns(), dataOptionDefaults()[["posCol"]])
  ),column(4,
           selectInput("accessionColType", "Accession Column Type", accTypes(), dataOptionDefaults()[["accColType"]]),
           selectInput("accessionCol", "Accession Column", dataColumns(), dataOptionDefaults()[["accCol"]]),
           selectInput("measurementCol", "Measurement Column", dataColumns(), dataOptionDefaults()[["measurementCol"]]),
           selectInput("replicateCol", "Replicate Column", dataColumns(), dataOptionDefaults()[["replicateCol"]]),
           selectInput("experimentCol", "Experiment Column", dataColumns(), dataOptionDefaults()[["expCol"]])
  ),column(4,
           checkboxInput("hasControls", "Are controls included?", TRUE),
           conditionalPanel(condition = "input.hasControls",
           selectInput("controlCol", "Control Column", dataColumns(), dataOptionDefaults()[["ctrlCol"]])),        
           uiOutput("uiOutput_controls")
  )
  )
  do.call(fluidRow, elements)
})

output$uiOutput_controls <- renderUI({    
  ctrlSelects <- list(
    selectInput("posCtrl", "Positive Control", ctrlTypes(), dataOptionDefaults()[["posCtrls"]], multiple=T),
    selectInput("negCtrl", "Negative Control", ctrlTypes(), dataOptionDefaults()[["negCtrls"]], multiple=T)
  )
  do.call(wellPanel, ctrlSelects)
})