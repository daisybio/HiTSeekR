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
  
  #if we are dealing with miRNAs add miRNA family name and ID
  if(input$screenType == "miRNA" && input$accessionColType == "MIMAT")
  {
    library(mirbase.db)
    fam <- as.data.frame(mirbaseFAMILY)
    mimat <- as.data.frame(mirbaseMATURE)
    result <- left_join(mimat, fam, by="mirna_id")
    data$Accession <- as.character(data$Accession)
    data <- left_join(data, result, by=c("Accession" = "mature_acc"))
    
    mergeRows <- function(y){
      if(length(unique(y)) > 1) return(paste(y, collapse="/"))
      else return(y[1])
    }
    
    data <- data %>% group_by(Experiment, Plate, Well.position) %>% summarise_each(funs(mergeRows))
    return(as.data.frame(data))
  }
  
  return(data)
})
