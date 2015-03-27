library(RCurl)
library(rjson)
library(foreach)

base64EncFile <- function(fileName){  
  return(base64(readChar(fileName, file.info(fileName)$size)))
}

setup.KPM <- function(list.of.indicator.matrices, 
                      algorithm="Greedy", strategy="GLONE", graphID=1,
                      removeBENs=FALSE, range, 
                      Kmin=0, Lmin=0, Kmax=0, Lmax=0, Kstep=1, Lstep=1, 
                      ATTACHED_TO_ID, 
                      computed.pathways=20, 
                      with.perturbation=FALSE){  
  
  #create a run id
  RNAice_RUN_ID <- paste(sample(c(LETTERS[1:6],0:9),6,replace=TRUE),collapse="")
  
  #base64 encode files
  datasetList <- datasetList.KPM(list.of.indicator.matrices,ATTACHED_TO_ID, RNAice_RUN_ID)
  
  # setup the json settings:
  KPMsettings <- toJSON(
    list(
      parameters=c(
        name=paste("RNAice run", RNAice_RUN_ID),
        algorithm=algorithm,
        strategy=strategy,
        removeBENs=tolower(as.character(removeBENs)),
        unmapped_nodes="Add to negative list",
        computed_pathways=computed.pathways,
        graphID=graphID,
        l_samePercentage="false",
        samePercentage_val=0,
        k_values=list(c(val=Kmin, val_step=Kstep, val_max=Kmax, use_range=tolower(as.character(range)), isPercentage="false")),
        l_values=list(
          c(val=Lmin, val_step=Lstep, val_max=Lmax, use_range=tolower(as.character(range)), isPercentage="false", datasetName=paste(RNAice_RUN_ID, 1, sep=""))
        )
        ), 
      withPerturbation=tolower(as.character(with.perturbation)),
      perturbation=list(c( # perturbation can be left out, if withPeturbation parameter is set to false.
        technique="Node-swap",
        startPercent=5,
        stepPercent=1,
        maxPercent=15,
        graphsPerStep=1
      )),      
      linkType="OR",
      attachedToID=ATTACHED_TO_ID,
      positiveNodes="",
      negativeNodes=""
    ))  
  return(list(KPMsettings, datasetList))
}

datasetList.KPM <- function(list.of.indicator.matrices, ATTACHED_TO_ID, RNAice_RUN_ID)
{  
  counter <- 0
  datasetList <- foreach(indicator.matrix = list.of.indicator.matrices) %do% {
    txt.con <- textConnection("tmp.file", "w")    
    
    write.table(indicator.matrix, txt.con, sep="\t",quote=F)    
    enc.file <- base64(tmp.file)
    close(txt.con)
    counter <- counter + 1
    c(name=paste(RNAice_RUN_ID, counter, sep=""), attachedToID=ATTACHED_TO_ID, contentBase64=enc.file)
  }
  
  return(toJSON(datasetList))
}

call.KPM <- function(indicator.matrices, ATTACHED_TO_ID=NULL, url="http://localhost:8080/kpm-web/",...){  
  # generate random UUID:
  if(is.null(ATTACHED_TO_ID))
  ATTACHED_TO_ID = paste(sample(c(LETTERS[1:6],0:9),32,replace=TRUE),collapse="")
  
  #PPI network for KPM if we wanted to use a non-default one
  #graphFile <- "data/biogrid_entrez.sif"
  #graphFile <- "data/graph-ulitsky-entrez.sif"
  #graphFile <- "data/graph-hprd-entrez.sif"
  #graph <- base64EncFile(graphFile)
  #graph <- toJSON(c(name=basename(graphFile), attachedToID=ATTACHED_TO_ID, contentBase64=graph))
  
  #KPM settings:
  kpmSetup <- setup.KPM(indicator.matrices, ATTACHED_TO_ID=ATTACHED_TO_ID, ...)
    
  result <- NULL
  print(sprintf("url: %s", url))
  print(sprintf("kpmSettings: %s", kpmSetup[[1]]))  
  result <- sendToKpmServiceAsync(url, kpmSetup, NULL)

  return(result)
}

withTryCatch <- function(surroundedFunc){
  tryCatch({
    surroundedFunc()
  }, error = function(e) {
    if("COULDNT_CONNECT" %in% class(e)){
      stop("Couldn't connect to KPM url.")
    }else{      
      stop(paste("Unexpected error in KPM:", class(e)))
    }    
    return(NULL)
  })
}

sendToKpmServiceAsync <- function(url, kpmSetup, inputGraph){
  withTryCatch(function(){    
    url <- paste(url, "requests/kpmAsyncJSON", sep="")    
    result <- postForm(url, kpmSettings=kpmSetup[[1]], datasets=kpmSetup[[2]])#, graph=inputGraph)
    
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
