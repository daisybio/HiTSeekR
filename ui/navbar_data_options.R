#identifier column types
accTypes <- reactive({
  if(is.null(input$screenType)) return(NULL)
  else if(input$screenType == "siRNA") return(c("RefSeq", "Entrez", "GeneSymbol", "FlybaseCG"))
  else if(input$screenType == "miRNA") return(c("MI", "MIMAT", "mature_name"))
  #return(c("RefSeq", "Entrez", "MIMAT", "MI"))
})

#Are controls included?
hasCtrls <- reactive({
  if(is.null(input$hasCtrls)) return(FALSE)
  else return(input$hasCtrls)
})

#return available control types found in data
ctrlTypes <- reactive({
  data <- rawData()
  if(!is.null(input$controlCol) && input$controlCol %in% colnames(data)) ctrlCol <- data[,input$controlCol]
  
  else ctrlCol <- data[,dataOptionDefaults()[["ctrlCol"]]]
  
  unique(as.character(ctrlCol))
})

#select a negative control
negCtrl <- reactive({
  if(is.null(input$negCtrl)){
    return(dataOptionDefaults()[["negCtrls"]])
  }
  else{
    input$negCtrl
  }
})

#select a positive control
posCtrl <- reactive({
  if(is.null(input$posCtrl)){
    return(dataOptionDefaults()[["posCtrls"]])
  }
  else{
    input$posCtrl
  }
})

#select corresponding columns in the dataset to process to a common format
output$uiOutput_data_options <- renderUI({
  
  elements <- list(column(4,
                          selectInput("screenType", "Type of screen", c("Gene (e.g. siRNA)" = "siRNA", "miRNA (e.g. inhibitor)" = "miRNA", "drugs", "Small compounds (e.g. drugs)")),
                          selectInput("sampleCol", "Sample Name Column", dataColumns()),
                          selectInput("plateCol", "Plate Column", dataColumns()),
                          selectInput("positionColType", "Position Column Type", c("Alpha well names" = "alpha", "Numeric" = "numeric", "Row / Column" = "rowcol")),
                          conditionalPanel("input.positionColType == 'rowcol'", 
                            selectInput("rowCol", "Row Column", dataColumns()),
                            selectInput("colCol", "Column Column", dataColumns())
                          ),
                          conditionalPanel("input.positionColType != 'rowcol'",
                            selectInput("positionCol", "Position Column", dataColumns())
                          )
  ),column(4,
           selectInput("accessionColType", "Accession Column Type", NULL),
           selectInput("accessionCol", "Accession Column", dataColumns()),           
           selectInput("measurementCol", "Measurement Column", dataColumns(), multiple=T),
           selectInput("replicateCol", "Replicate Column", dataColumns(), multiple=T),
           selectInput("experimentCol", "Experiment Column", dataColumns(), multiple=T)
  ),column(4,
           checkboxInput("hasControls", "Are controls included?", FALSE),
           conditionalPanel(condition = "input.hasControls",
           selectInput("controlCol", "Control Column", dataColumns()),        
           uiOutput("uiOutput_controls")),
           checkboxInput("log2normalize", "Log2 transform signal data", FALSE),
           checkboxInput("computeBscore", "Compute B-score", FALSE)
  )
  )
  do.call(fluidRow, elements)
})

#make an update of the data import options if necessary
dataUpdateObserver <- observe({
  datasetName()  
  isolate({
    updateSelectInput(session, "screenType", "Type of screen", c("Gene silencing" = "siRNA", "miRNA inhibitor / mimics" = "miRNA", "compound screen" = "compound"), dataOptionDefaults()[["screenType"]])
    updateSelectInput(session, "sampleCol", "Sample Name Column", dataColumns(), dataOptionDefaults()[["sampleCol"]])
    updateSelectInput(session, "plateCol", "Plate Column", dataColumns(), dataOptionDefaults()[["plateCol"]])
    updateSelectInput(session, "positionColType", "Position Column Type", c("Alpha well names" = "alpha", "Numeric" = "numeric", "Row / Column" = "rowcol"), dataOptionDefaults()[["posColType"]])
    updateSelectInput(session, "rowCol", "Row Column", dataColumns(), dataOptionDefaults()[["rowCol"]])
    updateSelectInput(session, "colCol", "Column Column", dataColumns(), dataOptionDefaults()[["colCol"]])
    updateSelectInput(session, "positionCol", "Position Column", dataColumns(), dataOptionDefaults()[["posCol"]])
    #updateSelectInput(session, "accessionColType", "Accession Column Type", accTypes(), dataOptionDefaults()[["accColType"]])
    updateSelectInput(session, "accessionCol", "Accession Column", dataColumns(), dataOptionDefaults()[["accCol"]])
    updateSelectInput(session, "measurementCol", "Measurement Column", dataColumns(), dataOptionDefaults()[["measurementCol"]])
    updateSelectInput(session, "replicateCol", "Replicate Column", dataColumns(), dataOptionDefaults()[["replicateCol"]])
    updateSelectInput(session, "experimentCol", "Experiment Column", dataColumns(), dataOptionDefaults()[["expCol"]])
    hasCtrls <- dataOptionDefaults()[["hasCtrls"]]=="TRUE"
    updateCheckboxInput(session, "hasControls", "Are controls included?", hasCtrls)
    if(hasCtrls){
      updateSelectInput(session, "controlCol", "Control Column", dataColumns(), dataOptionDefaults()[["ctrlCol"]])
      updateSelectInput(session, "posCtrl", "Positive Control", ctrlTypes(), dataOptionDefaults()[["posCtrls"]])
      updateSelectInput(session, "negCtrl", "Negative Control", ctrlTypes(), dataOptionDefaults()[["negCtrls"]])
    }
  })
}, priority=1000)


dataChangeObserver <- observe({
    datasetName()
    input$screenType
    input$sampleCol
    input$plateCol
    input$positionColType
    input$rowCol
    input$colCol
    input$positionCol
    input$accessionCol
    input$accessionColType
    input$measurementCol
    input$replicateCol
    input$experimentCol
    input$hasControls
    input$controlCol
    input$posCtrl
    input$negCtrl
    isolate({
      if(input$startButton != 0)
      showshinyalert(session, "general_status", "Input has changed. You need to process the data again for the changes to take effect", "danger")  
    })
})

accessionColTypeObserver <- observe({
  updateSelectInput(session, "accessionColType", "Accession Column Type", accTypes(), dataOptionDefaults()[["accColType"]])
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
#outputOptions(output, 'uiOutput_controls', suspendWhenHidden=FALSE)
