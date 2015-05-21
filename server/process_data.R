processedData <- reactive({  
  input$startButton

  isolate({
  hideshinyalert(session, "general_status")
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
  
  if(is.null(input$positionColType)) return(NULL)
  
  #first check if we are dealing with several readouts
  numberOfReadouts <- length(input$measurementCol)
  numberOfExperiments <- length(input$experimentCol)
  numberOfReplicates <- length(input$replicateCol)
  
  #count which types of data have multiple columns
  count <- 0
  if(numberOfReadouts > 1) count <- count + 1
  if(numberOfExperiments > 1) count <- count + 1
  if(numberOfReplicates > 1) count <- count + 1
  
  if(count > 1) stop("We do not support several experiments and measurements (readouts / replicates) in the same data set")
  
  if(numberOfReadouts == 0 && numberOfExperiments < 2 && numberOfReplicates < 2) stop("You need to specify at least one measurement column (e.g. fluorescence counts)")
  
  if(numberOfReadouts  == 1){
    data$readoutCol <- input$measurementCol
    data$measurementCol  <- data[,input$measurementCol]  
  }
  else{
    data <- gather_(data, "readoutCol", "measurementCol", input$measurementCol) 
  }
    
  if(numberOfExperiments == 0) data$experimentCol <- datasetName()
  else if(numberOfExperiments == 1) data$experimentCol <- data[,input$experimentCol]
  else{
    data <- gather_(data, "experimentCol", "measurementCol", input$experimentCol)
  }
  
  if(numberOfReplicates == 0) data$replicateCol <- 1
  else if(numberOfReplicates == 1) data$replicateCol <- data[,input$replicateCol]
  else{
    data <- gather_(data, "replicateCol", "measurementCol", input$replicateCol)
  }
  
  sampleCol <- data[,input$sampleCol]  
  
  if(input$positionColType == "alpha")
  {
    wellAlpha <- repairAlphaName(data[,input$positionCol])
    rowCol <- alphaNames2Pos(wellAlpha)
  }
  else if(input$positionColType == "numeric")
  {
    wellAlpha <- positionsToAlphaName(as.integer(data[,input$positionCol]))
    rowCol <- alphaNames2Pos(wellAlpha)
  }
  else if(input$positionColType == "rowcol")
  {
    rows <- data[,input$rowCol]
    cols <- data[,input$colCol]
    rowCol <- cbind(rows, cols)
    wellAlpha <- foreach(row = iter(data, "row"), .inorder = TRUE, .combine=c) %do% {
      paste(LETTERS[row[,input$colCol]], row[,input$rowCol], sep="")
    }
  }
  accessionCol <- data[,input$accessionCol]
  plateCol <- data[,input$plateCol]
  
  if(input$hasControls) controlCol <- data[,input$controlCol]
  else controlCol <- NA
  
  processedData <- data.frame(data$experimentCol, sampleCol, accessionCol, plateCol, wellAlpha, rowCol, data$replicateCol, controlCol, data$readoutCol, data$measurementCol)
  colnames(processedData) <- c("Experiment", "Sample", "Accession", "Plate", "Well.position", "Row", "Column", "Replicate", "Control", "Readout", "Raw")
  
  if(input$log2normalize){
    processedData$Raw <- log2(processedData$Raw)
  }
  
  if(input$hasControls)
  {  
    result <- normalizeRawData(processedData, control.based=T, pos.ctrl=input$posCtrl, neg.ctrl=input$negCtrl, updateProgress=updateProgress, compute.B=input$computeBscore)
  }
  else{
    result <- normalizeRawData(processedData, control.based=F, updateProgress=updateProgress, compute.B=input$computeBscore)    
  }
  return(as.data.frame(result))
  })
})
