#filter and summarise
data <- reactive({
  
  data <- processedData()
  
  progress <- shiny::Progress$new()
  progress$set(message = "Merge & Filter...", value = 0)
  on.exit(progress$close())
  
  #update progress bar function
  updateProgress <- function(value = NULL, detail = NULL) {
    if (is.null(value)) {
      value <- progress$getValue()
      value <- value + (progress$getMax() - value) / 5
    }
    progress$set(value = value, detail = detail)
  }
  
  input$updateExclusion
  
  data <- isolate({
    if(input$updateExclusion != 0 && nchar(input$exclude) > 0){
      data <- data[-grep(input$exclude, data$Sample),]
    }
    else data
  })
  
  updateProgress(detail = "Merging replicates", value=0.2)
  deviations <- ddply(data, .(Experiment, Plate, Well.position), numcolwise(function(x){sd(x)/length(x)}))
  data <- ddply(data, .(Experiment, Plate, Well.position,Sample, Accession, Control), numcolwise(mean))
  data <- merge(data, deviations[,setdiff(colnames(deviations),c("wellCount", "Row", "Column", "Sample", "Control"))], by=c("Experiment", "Plate", "Well.position"), suffixes=c("", ".sem"))
  
  #if we are dealing with miRNAs add miRNA family name and ID
  if(input$screenType == "miRNA")
  {
    updateProgress(detail = "Querying mirbase", value=0.6)
    library(mirbase.db)
    fam <- as.data.frame(mirbaseFAMILY)
    
    mergeRows <- function(y){
      if(length(unique(y)) > 1) return(paste(y, collapse="/"))
      else return(y[1])
    }     
    data$Accession <- as.character(data$Accession)
    
    if(input$accessionColType == "MIMAT")
    {
      mimat <- as.data.frame(mirbaseMATURE)
      result <- left_join(mimat, fam, by="mirna_id")
      data <- left_join(data, result, by=c("Accession" = "mature_acc"))
      data <- data %>% group_by(Experiment, Plate, Well.position) %>% summarise_each(funs(mergeRows))      
    }
    
    else if(input$accessionColType=="MI"){
      mi <- as.data.frame(mirbaseACC2ID)
      data <- left_join(data, mi, by=c("Accession" = "mirna_acc"))
      data <- left_join(data, fam, by=c("mirna_id"))
    }
    
    else{
      stop("unknown miRNA identifier selected")
    }
    
    return(as.data.frame(data))
  }
  
  return(data)
})
