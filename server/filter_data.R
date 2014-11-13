#filter and summarise
data <- reactive({
  
  data <- processedData()
  
  signalType <- input$normalization
  data$signal <- data[[signalType]]
  input$updateExclusion
  
  data <- isolate({
    if(input$updateExclusion != 0 && nchar(input$exclude) > 0){
      data <- data[-grep(input$exclude, data$miRBase.ID.miRPlus.ID),]
    }
    else data
  })
  deviations <- ddply(data, .(Plate, Well.position), numcolwise(function(x){sd(x)/length(x)}))
  data <- ddply(data, .(Plate, Well.position,miRBase.ID.miRPlus.ID, miRBase.accession, Sample, Control), numcolwise(mean))
  data <- merge(data, deviations[,setdiff(colnames(deviations),c("wellCount", "Row", "Column", "Sample", "Control"))], by=c("Plate", "Well.position"), suffixes=c("", ".sem"))
  return(data)
})
