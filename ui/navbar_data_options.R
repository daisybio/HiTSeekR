#identifier column types
accTypes <- reactive({
  if(is.null(input$screenType)) return(NULL)
  else if(input$screenType == "siRNA") return(c("RefSeq", "Entrez", "GeneSymbol"))
  else if(input$screenType == "miRNA") return(c("MI", "MIMAT"))
  #return(c("RefSeq", "Entrez", "MIMAT", "MI"))
})

#return available control types found in data
ctrlTypes <- reactive({
  data <- rawData()
  if(!is.null(input$controlCol)) ctrlCol <- data[,input$controlCol]
  
  else ctrlCol <- data[,dataOptionDefaults()[["controlCol"]]]
  
  unique(as.character(ctrlCol))
})

#select a negative control
negCtrl <- reactive({
  if(is.null(input$negCtrl)){
    return(dataOptionDefaults()[["negCtrls"]])
  }
})

#select a positive control
posCtrl <- reactive({
  if(is.null(input$posCtrl)){
    return(dataOptionDefaults()[["posCtrls"]])
  }
})

#select corresponding columns in the dataset to process to a common format
output$uiOutput_data_options <- renderUI({
  elements <- list(column(4,
                          selectInput("screenType", "Type of screen", c("Gene knockout (e.g. siRNA)" = "siRNA", "miRNA inhibitor / mimics" = "miRNA"), dataOptionDefaults()[["screenType"]]),
                          selectInput("sampleCol", "Sample Name Column", dataColumns(), dataOptionDefaults()[["sampleCol"]]),
                          selectInput("plateCol", "Plate Column", dataColumns(), dataOptionDefaults()[["plateCol"]]),
                          selectInput("positionColType", "Position Column Type", c("Alpha well names" = "alpha", "Numeric" = "numeric", "Row / Column" = "rowcol"), dataOptionDefaults()[["posColType"]]),
                          conditionalPanel("input.positionColType == 'rowcol'", 
                            selectInput("rowCol", "Row Column", dataColumns(), dataOptionDefaults()[["rowCol"]]),
                            selectInput("colCol", "Column Column", dataColumns(), dataOptionDefaults()[["colCol"]])
                          ),
                          conditionalPanel("input.positionColType != 'rowcol'",
                            selectInput("positionCol", "Position Column", dataColumns(), dataOptionDefaults()[["posCol"]])
                          )
  ),column(4,
           selectInput("accessionColType", "Accession Column Type", accTypes(), dataOptionDefaults()[["accColType"]]),
           selectInput("accessionCol", "Accession Column", dataColumns(), dataOptionDefaults()[["accCol"]]),
           selectInput("measurementCol", "Measurement Column", dataColumns(), dataOptionDefaults()[["measurementCol"]]),
           selectInput("replicateCol", "Replicate Column", dataColumns(), dataOptionDefaults()[["replicateCol"]]),
           selectInput("experimentCol", "Experiment Column", dataColumns(), dataOptionDefaults()[["expCol"]])
  ),column(4,
           checkboxInput("hasControls", "Are controls included?", dataOptionDefaults()[["hasCtrls"]]=="TRUE"),
           conditionalPanel(condition = "input.hasControls",
           selectInput("controlCol", "Control Column", dataColumns(), dataOptionDefaults()[["ctrlCol"]]),        
           uiOutput("uiOutput_controls"))
  )
  )
  do.call(fluidRow, elements)
})

#extra well pane for selecting controls of the experiment
output$uiOutput_controls <- renderUI({    
  ctrlSelects <- list(
    selectInput("posCtrl", "Positive Control", ctrlTypes(), dataOptionDefaults()[["posCtrls"]], multiple=T),
    selectInput("negCtrl", "Negative Control", ctrlTypes(), dataOptionDefaults()[["negCtrls"]], multiple=T)
  )
  do.call(wellPanel, ctrlSelects)
})

#This bit is important to make sure that the file import settings are updated
#with the defaults of the selected dataset even though the corresponding input
#is hidden in the UI.
outputOptions(output, 'uiOutput_data_options', suspendWhenHidden=FALSE)
