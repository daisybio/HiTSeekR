#load data
rawData <- reactive({
  
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  
  if(input$getFromPubChemButton > 0)
  {    
    progress$set(message = "Downloading data from PubChem...", value=0)
    tryCatch({
      library(RCurl)
      raw <- getURL(paste("https://pubchem.ncbi.nlm.nih.gov/assay/assay.cgi?aid=", input$pubchem_aid, "&q=expdata_csvsave", sep=""))
      if(grepl("err message:", raw)){
        stop(raw)
      }
      data <- read.csv(text = raw)
    }, 
    warning = function(w){ showshinyalert(session, "general_status", paste("potential problem when reading file:", w$message), "warning")},
    error = function(e){ showshinyalert(session, "general_status", paste("error reading file:", e$message), "danger") })
  }
  if(is.null(input$file)){
    if(input$dataset != "none selected") progress$set(message = paste("Loading demo data set", input$dataset), value=0)
    if(input$dataset == "BCSC")
      data <- read.delim(paste(data.folder, "BCSC.txt", sep=""), header=T)
    else if(input$dataset == "A375_MTS"){
      data <- read.table(paste(data.folder, "A375_MTS.txt",sep=""), header=T, sep="\t")
    }
    else if(input$dataset == "HCC_vorinostat_siRNA")
    {
      data <- read.csv(paste(data.folder, "AID_743454_data.csv", sep=""), header=T)
    }
    else if(input$dataset == "HCC_vorinostat_miRNA")
    {
      data <- read.csv(paste(data.folder, "AID_743456_MIMAT_data.csv", sep=""), header=T)
    }
    else if(input$dataset == "TNFa_Casp4")
    {
      data <- read.delim(paste(data.folder, "Caspase4_processed.txt", sep=""), header=T, sep="\t")            
    }
    else if(input$dataset == "KRAS_synleth_compound"){
      data <- read.delim(paste(data.folder, "KRasSyntheticLethal_CellTiterGlo(1054.0014)_screeningdata.txt", sep=""), header=T, sep="\t", fill=T)
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
    progress$set(message = "Uploading and processing custom data set...", value=0)
    file <- input$file
    tryCatch({
    data <- read.table(file$datapath, header=T, sep=input$fileSeparator, fill=T)
    }, 
    warning = function(w){ showshinyalert(session, "general_status", paste("potential problem when reading file:", w$message), "warning")},
    error = function(e){ showshinyalert(session, "general_status", paste("error reading file:", e$message), "danger") })
  }
  progress$set(message =  "Loading complete", value=1)
  return(data)
})

datasetName <- reactive({
  if(!is.null(input$file) || input$getFromPubChemButton > 0){
    updateSelectInput(session, "dataset", c(demo.data.sets, "CUSTOM"="CUSTOM"), "CUSTOM")
    return("CUSTOM")
  } 
  else return(input$dataset)
})
