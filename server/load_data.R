#load data
rawData <- reactive({
  
  if(is.null(input$file)){
    if(input$dataset == "BCSC")
      data <- read.delim("data/BCSC.txt", header=T)
    else if(input$dataset == "A375_MTS"){
      data <- read.table("data/A375_MTS.txt", header=T, sep="\t")
    }
  }
  else{
    file <- input$file
    data <- read.table(file$datapath, header=T, sep=input$fileSeparator)
  }
  
  return(data)
})
