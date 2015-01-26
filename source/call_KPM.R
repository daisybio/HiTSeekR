library(RCurl)
library(rjson)
library(foreach)

base64EncFile <- function(fileName){  
  return(base64(readChar(fileName, file.info(fileName)$size)))
}

setup.KPM <- function(list.of.indicator.matrices, graphFile, algorithm="Greedy", strategy="GLONE", removeBENs=FALSE, K=0, L=0, ATTACHED_TO_ID){
  
  #create tempfiles
  tmp.file <- tempfile()
  
  #base64 encode files
  datasetList <- datasetList.KPM(list.of.indicator.matrices, list(tmp.file), ATTACHED_TO_ID)
  
  #create a run id
  RNAice_RUN_ID <- paste(sample(c(LETTERS[1:6],0:9),6,replace=TRUE),collapse="")
  
  # setup the json settings:
  KPMsettings <- toJSON(
    list(
      parameters=c(
        name=paste("RNAice run", RNAice_RUN_ID),
        algorithm=algorithm,
        strategy=strategy,
        removeBENs=tolower(as.character(removeBENs)),
        unmapped_nodes="Add to positive list",
        computed_pathways=20,
        graphName=basename(graphFile),             
        l_samePercentage="false",
        samePercentage_val=0,
        k_values=list(c(val=K, val_step=1, val_max=1, use_range="false", isPercentage="false")),
        l_values=list(      
          c(val=L, val_step=1, val_max=1, use_range="false", isPercentage="false", datasetName=basename(tmp.file))
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

call.KPM <- function(indicator.matrices, url=url <- "http://localhost:8080/kpm-web/", ATTACHED_TO_ID=NULL, ...){
  
  # generate random UUID:
  if(is.null(ATTACHED_TO_ID))
  ATTACHED_TO_ID = paste(sample(c(LETTERS[1:6],0:9),32,replace=TRUE),collapse="")
  print(ATTACHED_TO_ID)
  
  #PPI network for KPM
  #graphFile <- "data/biogrid_entrez.sif"
  #graphFile <- "data/graph-ulitsky-entrez.sif"
  graphFile <- "data/graph-hprd-entrez.sif"
  graph <- base64EncFile(graphFile)
  graph <- toJSON(c(name=basename(graphFile), attachedToID=ATTACHED_TO_ID, contentBase64=graph))
  
  #KPM settings:
  kpmSetup <- setup.KPM(indicator.matrices, ATTACHED_TO_ID=ATTACHED_TO_ID, graphFile, ...)
    
  result <- NULL
  print(sprintf("url: %s", url))
  print(sprintf("kpmSettings: %s", kpmSetup[[1]]))  
  result <- sendToKpmServiceAsync(url, kpmSetup, graph)

  return(result)
}

withTryCatch <- function(surroundedFunc){
  tryCatch({
    surroundedFunc()
  }, error = function(e) {
    if("COULDNT_CONNECT" %in% class(e)){
      print("Error: Couldn't connect to url.")
    }else{
      print("Unexpected error:")
      print(class(e))
    }    
    return(NULL)
  })
}

sendToKpmServiceAsync <- function(url, kpmSetup, inputGraph){
  withTryCatch(function(){    
    url <- paste(url, "requests/kpmAsyncJSON", sep="")    
    result <- postForm(url, kpmSettings=kpmSetup[[1]], datasets=kpmSetup[[2]], graph=inputGraph)
    
    jsonResult <- fromJSON(result)
    print(jsonResult["comment"])   
    return(jsonResult)
  })
}

getKpmRunStatus <- function(url, questId){
  withTryCatch(function(){
    url <- paste(url, "requests/kpmRunStatus", sep="")
    print(sprintf("url: %s", url))    
    result <- postForm(url, questID=questId)
    jsonResult <- fromJSON(result)
    
    if(tolower(jsonResult["success"]) == "cancelled"){
      print("Run has been cancelled.")
      return
    }
    
    print(jsonResult["completed"]) 
    print(jsonResult["progress"])
    
    return(jsonResult)
  })
}

getKpmResults <- function(url, questId){
  withTryCatch(function(){
    url <- paste(url, "requests/kpmResults", sep="")
    print(sprintf("url: %s", url))
    
    result <- postForm(url, questID=questId)
    jsonResult <- fromJSON(result)
    
    if(tolower(jsonResult["success"]) == "true"){
      return(jsonResult)
    }
    else{
      return(NULL)
    }
  })
}
