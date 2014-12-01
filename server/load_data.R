#load data
rawData <- reactive({
  
  if(is.null(input$file)){
    if(input$dataset == "BCSC")
      data <- read.delim("data/BCSC.txt", header=T)
    else if(input$dataset == "MTS data"){
      data <- read.table("data/jos_processed.txt", header=T, sep="\t")
      #data$Control <- NA
      #data$Sample <- NA
      #data$Accession <- NA
    }
  }
  else{
    file <- input$file
    data <- read.table(file$datapath, header=T, sep=input$fileSeparator)
  }
  
  #data$Plate <- as.factor(data$Plate)
  #data$Replicate <- as.factor(data$Replicate)
  #shorten miRNA label names
  #data$Sample <- gsub(",.*", "", gsub("hsa-", "", data$Sample))
  
  return(data)
})
