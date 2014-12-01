#filter and summarise
data <- reactive({
  
  data <- processedData()
  
  input$updateExclusion
  
  data <- isolate({
    if(input$updateExclusion != 0 && nchar(input$exclude) > 0){
      data <- data[-grep(input$exclude, data$Sample),]
    }
    else data
  })
  deviations <- ddply(data, .(Experiment, Plate, Well.position), numcolwise(function(x){sd(x)/length(x)}))
  data <- ddply(data, .(Experiment, Plate, Well.position,Sample, Accession, Control), numcolwise(mean))
  data <- merge(data, deviations[,setdiff(colnames(deviations),c("wellCount", "Row", "Column", "Sample", "Control"))], by=c("Experiment", "Plate", "Well.position"), suffixes=c("", ".sem"))
  return(data)
})
