#per session attached to id for KPM
ATTACHED_TO_ID <- paste(sample(c(LETTERS[1:6],0:9),32,replace=TRUE),collapse="")

KPM.run <- reactive({
  input$startKPMButton
  isolate({
    showshinyalert(session, "kpm_status", "Generating indicator matrix", "info")  
    indicator.matrix <- targets.indicator.matrix()
    showshinyalert(session, "kpm_status", "Sending data to KPM", "info")
    #browser()
    if(input$kpm_ranged)
    {
      Kmin <- input$kpm_lower_K
      Lmin <- input$kpm_lower_L
    }
    else{
      Kmin <- input$kpm_K
      Lmin <- input$kpm_L
    }
    result <- call.KPM(list(indicator.matrix), 
                       url=input$kpm_URL, 
                       ATTACHED_TO_ID=ATTACHED_TO_ID, 
                       strategy=input$kpm_strategy, 
                       algorithm=input$kpm_algorithm, 
                       removeBENs=input$kpm_ben_removal, 
                       range=input$kpm_ranged,
                       Kmin=Kmin, 
                       Lmin=Lmin,
                       Kmax=input$kpm_upper_L,
                       Lmax=input$kpm_upper_L,
                       Kstep=input$kpm_step_K,
                       Lstep=input$kpm_step_L,
                       with.perturbation=input$kpm_perturbation,
                       computed.pathways=input$kpm_pathways)
    return(result)
  })
})

output$KPM.test <- renderPrint({
  
  if(input$startKPMButton == 0){
    showshinyalert(session, "kpm_status", "Press the start button to initiate a remote KeyPathwayMiner analysis", "info")      
    return(NULL)
  }
  else{
    result <- KPM.result()
  }
  
  return(print(result))
})

quest.progress.url <- function(){
  kpm.url <- paste(input$kpm_URL, "requests/quests?attachedToId=", sep="")
  paste("Click <a target='_blank' href='", kpm.url, ATTACHED_TO_ID, "&hideTitle=false", "'>here</a> to follow the progress of your run in KPM web", sep="")
}

quest.result.url <- function(){
  kpm.url <- paste(input$kpm_URL, "requests/quests?attachedToId=", sep="")
  paste("Click <a target='_blank' href='", kpm.url, ATTACHED_TO_ID, "&hideTitle=false", "'>here</a> to see the results of your run in KPM web", sep="")  
}

check.KPM.result <- function(){

  #KPM hasn't been executed yet
  if(input$startKPMButton == 0) return(FALSE)
  
  #get current running status from KPM web
  quest.status <- getKpmRunStatus(input$kpm_URL, currentQuest())
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
  return(getKpmResults(input$kpm_URL, currentQuest()))
}

currentQuest <- reactive({
  start.run <- KPM.run()
  return(start.run[["questID"]])
})

KPM.result <- reactive({
  KPM.finished <- check.KPM.result()
  if(KPM.finished) return(get.KPM.result())
  else{
    invalidateLater(5000, session)
    return(NULL)
  }
})