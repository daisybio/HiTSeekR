library(RCurl)
library(rjson)

base64EncFile <- function(fileName){  
  return(base64(readChar(fileName, file.info(fileName)$size)))
}

setup.KPM <- function(list.of.indicator.matrices, graphFile, algorithm="Greedy", strategy="INES", removeBENs=FALSE, ATTACHED_TO_ID){
  
  #create tempfiles
  tmp.file <- tempfile()
  
  #base64 encode files
  datasetList <- datasetList.KPM(list.of.indicator.matrices, list(tmp.file), ATTACHED_TO_ID)
  
  # setup the json settings:
  KPMsettings <- toJSON(
    list(
      parameters=c(
        name=paste("RNAice run", ATTACHED_TO_ID),
        algorithm=algorithm,
        strategy=strategy,
        removeBENs=as.character(removeBENs),
        unmapped_nodes="Add to positive list",
        computed_pathways=20,
        graphName=basename(graphFile),             
        l_samePercentage="false",
        samePercentage_val=0,
        k_values=list(c(val=5, val_step=1, val_max=1, use_range="false", isPercentage="false")),
        l_values=list(      
          c(val=102, val_step=1, val_max=1, use_range="false", isPercentage="false", datasetName=basename(tmp.file))
        )), 
      withPerturbation="false",
      linkType="OR",
      attachedToID=ATTACHED_TO_ID,
      positiveNodes="",
      negativeNodes=""
    ))
  return(list(KPMsettings, datasetList))
}

datasetList.KPM <- function(list.of.indicator.matrices, list.of.tmp.files, ATTACHED_TO_ID)
{
  counter <- 1
  datasetList <- foreach(indicator.matrix = list.of.indicator.matrices) %do% {
    tmp.file <- list.of.tmp.files[[counter]]
    counter <- counter + 1
    
    write.table(indicator.matrix, tmp.file, sep="\t",quote=F)    
    enc.file <- base64EncFile(tmp.file)
    
    c(name=basename(tmp.file), attachedToID=ATTACHED_TO_ID, contentBase64=enc.file)
  }
  
  return(toJSON(datasetList))
}

call.KPM <- function(indicator.matrices){
  
  # generate random UUID:
  ATTACHED_TO_ID = paste(sample(c(letters[1:6],0:9),30,replace=TRUE),collapse="")
  
  #PPI network for KPM
  #graphFile <- "data/biogrid_entrez.sif"
  graphFile <- "data/graph-ulitsky-entrez.sif"
  graph <- base64EncFile(graphFile)
  graph <- toJSON(c(name=basename(graphFile), attachedToID=ATTACHED_TO_ID, contentBase64=graph))
  
  #KPM settings:
  kpmSetup <- setup.KPM(indicator.matrices, ATTACHED_TO_ID=ATTACHED_TO_ID, graphFile)
    
  result <- NULL
  tryCatch({
    
    url <- "http://localhost:8080/kpm-web/requests/kpmJSON"
    print(sprintf("url: %s", url))
    
    result <- postForm(url, kpmSettings=kpmSetup[[1]], datasets=kpmSetup[[2]], graph=graph)
    result <- fromJSON(result)  
  }, error = function(e) {
    
    if("COULDNT_CONNECT" %in% class(e)){
      print("Error: Couldn't connect to url.")
    }else{
      print("Unexpected error:")
      print(class(e))
    }
    
  })  
  return(result)
}