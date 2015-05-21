#per session attached to id for KPM
ATTACHED_TO_ID <- paste(sample(c(LETTERS[1:6],0:9),32,replace=TRUE),collapse="")

KPM.network.list <- reactive({  
  tryCatch({
    kpm.url <- paste(keypathwayminer.url, "requests/graphsAsJSON/", sep="")    
    result <- getURL(kpm.url)
    jsonResult <- fromJSON(result)
    networks <- foreach(network = jsonResult, .combine=append) %do% {network[[1]]} 
    names(networks) <- foreach(network = jsonResult, .combine=append) %do% {network[[2]]} 
    return(networks)
  }, error = function(e) return(NULL))  
})

KPM.indicator.matrix <- reactive({
  if(input$screenType == "miRNA")              
  {                     
    return(targets.indicator.matrix())
  } else 
  {        
    return(genes.indicator.matrix())
  }     
})

KPM.run <- reactive({
  if(input$startKPMButton == 0) return(NULL)
  isolate({
    #if(input$kpm_ranged && input$random.miRNA.test) stop("miRNA target permutation test is limited to specific K and L")
    showshinyalert(session, "kpm_status", "Generating indicator matrix", "info")      
    indicator.matrix <- KPM.indicator.matrix()
    
    list.of.ind.matrices <- list(indicator.matrix)
    
    #if(input$random.miRNA.test){
    #  showshinyalert(session, "kpm_status", "Generating matrices for random miRNAs", "info")    
    #  list.of.ind.matrices <- append(list.of.ind.matrices, list.of.random.mirna.indicator.matrices())
    #}      
    
    showshinyalert(session, "kpm_status", "Sending data to KPM", "info")
    
    #if(input$kpm_ranged)
    #{
    #  Kmin <- input$kpm_lower_K
    #  Lmin <- input$kpm_lower_L
    #}
    #else{
      Kmin <- input$kpm_K
      Lmin <- input$kpm_L
    #}
    
    result <- foreach(ind.matrix = list.of.ind.matrices) %do%{
                    call.KPM(list(ind.matrix), 
                       url=keypathwayminer.url, 
                       ATTACHED_TO_ID=ATTACHED_TO_ID, 
                       strategy=input$kpm_strategy, 
                       algorithm=input$kpm_algorithm, 
                       graphID=input$kpm_network,
                       removeBENs=input$kpm_ben_removal, 
                       range=FALSE,
                       #range=input$kpm_ranged,
                       Kmin=Kmin, 
                       Lmin=Lmin,
                       #Kmax=input$kpm_upper_L,
                       #Lmax=input$kpm_upper_L,
                       #Kstep=input$kpm_step_K,
                       #Lstep=input$kpm_step_L,
                       #with.perturbation=input$kpm_perturbation,
                       computed.pathways=input$kpm_pathways)
              }
    #kpm.result <<- result
    return(result[[1]])
  })
})

output$KPM.test <- renderPrint({
  
  if(input$startKPMButton == 0){         
    return(NULL)
  }
  else{
    result <- KPM.result()
  }
  
  return(print(result))
})

quest.progress.url <- function(){
  kpm.url <- paste(keypathwayminer.url, "requests/quests?attachedToId=", sep="")
  paste("Click <a target='_blank' href='", kpm.url, ATTACHED_TO_ID, "&hideTitle=false", "'>here</a> to follow the progress of your run in KPM web", sep="")
}

quest.result.url <- function(){
  kpm.url <- paste(keypathwayminer.url, "requests/quests?attachedToId=", sep="")
  paste("Click <a target='_blank' href='", kpm.url, ATTACHED_TO_ID, "&hideTitle=false", "'>here</a> to see the results of your run in KPM web", sep="")  
}

check.KPM.result <- function(){

  #KPM hasn't been executed yet
  if(input$startKPMButton == 0){    
    return(FALSE)
  } 
  
  #get current running status from KPM web
  quest.status <- getKpmRunStatus(keypathwayminer.url, currentQuest())
  if(is.null(quest.status)) return(FALSE)
  
  #update status
  showshinyalert(session, "kpm_status", paste("KPM run is ", quest.status[["progress"]]*100, "% completed.", quest.progress.url()), "info")  
  
  if(quest.status[["completed"]] == TRUE){
    #update quest.completed to avoid unneccessary polling    
    showshinyalert(session, "kpm_status", paste("KPM run finished", quest.result.url()), "info")  
    
    return(quest.status[["completed"]])  
  }
  else{
    return(FALSE)
  }
}

get.KPM.result <- function(){
  return(getKpmResults(keypathwayminer.url, currentQuest()))
}

currentQuest <- reactive({
  start.run <- KPM.run()
  return(start.run[["questID"]])
})

KPM.result <- reactive({
  if(input$startKPMButton == 0){
    showshinyalert(session, "kpm_status", "Press the start button to initiate a remote KeyPathwayMiner analysis", "info") 
    return(NULL)
  } 
  KPM.finished <- check.KPM.result()
  if(KPM.finished) return(get.KPM.result())
  else{
    invalidateLater(5000, session)
    return(NULL)
  }
})