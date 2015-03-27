#load data
rawData <- reactive({
  if(input$getFromPubChemButton > 0)
  {
    library(RCurl)
    raw <- getURL(paste("https://pubchem.ncbi.nlm.nih.gov/assay/assay.cgi?aid=", input$pubchem_aid, "&q=expdata_csvsave", sep=""))
    data <- read.csv(text = raw)
  }
  if(is.null(input$file)){
    if(input$dataset == "BCSC")
      data <- read.delim("data/BCSC.txt", header=T)
    else if(input$dataset == "A375_MTS"){
      data <- read.table("data/A375_MTS.txt", header=T, sep="\t")
    }
    else if(input$dataset == "DM_Kc167"){
      library(cellHTS2)
      experimentName <- "KcViab"
      dataPath <- system.file(experimentName, package="cellHTS2")
      x <- readPlateList("Platelist.txt", name = experimentName,path = dataPath)
      x <- configure(x, descripFile = "Description.txt", confFile = "Plateconf.txt", logFile = "Screenlog.txt", path = dataPath)
      x <- cellHTS2::annotate(x, geneIDFile = "GeneIDs_Dm_HFA_1.1.txt", path = dataPath)
      kcViab <- cbind(as(featureData(x), "data.frame"), assayData(x)[["Channel 1"]])
      kcViab <- melt(kcViab,id.vars=c("plate", "well", "controlStatus", "HFAID", "GeneID"))
      colnames(kcViab)[ncol(kcViab)-1] <- "replicate"   
      data <- kcViab
      data$experiment <- experimentName
    }
  }
  else{
    file <- input$file
    data <- read.table(file$datapath, header=T, sep=input$fileSeparator)
  }
  
  return(data)
})

datasetName <- reactive({
  if(!is.null(input$file) || input$getFromPubChemButton > 0){
    updateSelectInput(session, "dataset", c(demo.data.sets, "CUSTOM"="CUSTOM"), "CUSTOM")
    return("CUSTOM")
  } 
  else return(input$dataset)
})
