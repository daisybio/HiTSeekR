processedData <- reactive({
  data <- rawData()
  
  progress <- shiny::Progress$new()
  progress$set(message = "Processing...", value = 0)
  on.exit(progress$close())
  
  #update progress bar function
  updateProgress <- function(value = NULL, detail = NULL) {
    if (is.null(value)) {
      value <- progress$getValue()
      value <- value + (progress$getMax() - value) / 5
    }
    progress$set(value = value, detail = detail)
  }
  
  sampleCol <- data[,input$sampleCol]
  if(input$positionColType == "alpha")
  {
    wellAlpha <- repairAlphaName(data[,input$positionCol])
    rowCol <- alphaNames2Pos(wellAlpha)
  }
  else if(input$positionColType == "numeric")
  {
    #TODO
  }
  replicateCol <- data[,input$replicateCol]
  measurementCol <- data[,input$measurementCol]
  accessionCol <- data[,input$accessionCol]
  plateCol <- data[,input$plateCol]
  experimentCol <- data[,input$experimentCol]
  controlCol <- data[,input$controlCol]
  
  processedData <- cbind(experimentCol, sampleCol, accessionCol, plateCol, wellAlpha, rowCol, replicateCol, controlCol, measurementCol)
  colnames(processedData) <- c("Experiment", "Sample", "Accession", "Plate", "Well.position", "Row", "Column", "Replicate", "Control", "Raw")
  
  if(input$hasControls && !(is.null(input$controlCol)))
  {  
    result <- normalizeRawData(processedData, control.based=T, pos.ctrl=input$posCtrl, neg.ctrl=input$negCtrl, updateProgress=updateProgress)
  }
  else if(!is.null(negCtrl()) && !is.null(posCtrl)){
    result <- normalizeRawData(processedData, control.based=T, pos.ctrl=posCtrl(), neg.ctrl=negCtrl(), updateProgress=updateProgress)    
  }
  else{
    result <- normalizeRawData(processedData, control.based=F, updateProgress=updateProgress)    
  }
  
  return(as.data.frame(result))
})
