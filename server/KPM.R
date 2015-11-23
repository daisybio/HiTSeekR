#per session attached to id for KPM
ATTACHED_TO_ID <- paste(sample(c(LETTERS[1:6],0:9),32,replace=TRUE),collapse="")

KPM.network.list <- reactive({  
  #create a progress bar
  progress <- shiny::Progress$new()
  progress$set(message = "Connecting to KeyPathwayMiner web to learn about available networks...", value = 0)
  on.exit(progress$close())
  networks <- NULL
  
  tryCatch({
    kpm.url <- paste(keypathwayminer.url, "rest/availableNetworks/", sep="")    
    result <- getURL(kpm.url)
    jsonResult <- fromJSON(result)
    networks <- foreach(network = jsonResult, .combine=append) %do% {network[[1]]} 
    names(networks) <- foreach(network = jsonResult, .combine=append) %do% {network[[2]]} 
    return(networks)
  }, error = function(e) {
      showshinyalert(session, "kpm_status", "Could not connect to KeyPathwayMiner Web", "danger")   
    })
  return(networks)
})

KPM.indicator.matrix <- reactive({
  if(input$screenType == "miRNA")              
  {                     
    return(targets.indicator.matrix())
  } else if(input$screenType == "compound")
  {
    return(drug.indicator.matrix())
  } else
  {        
    return(genes.indicator.matrix())
  }     
})

genes.indicator.matrix <- reactive({
  hits <- outliers()
  all.samples <- data()
  
  gene_ids <- na.omit(unique(all.samples$gene_id))
  gene_ids <- gene_ids[which(gene_ids %in% hits$gene_id)]
  ind.matrix <- as.matrix(rep(1, length(gene_ids)))
  row.names(ind.matrix) <- gene_ids
  
  return(ind.matrix)
})

KPM.run <- reactive({
  if(input$startKPMButton == 0) return(NULL)
  isolate({
    #if(input$kpm_ranged && input$random.miRNA.test) stop("miRNA target permutation test is limited to specific K and L")
    showshinyalert(session, "kpm_status", "Generating indicator matrix", "info")   
    indicator.matrix <- KPM.indicator.matrix()
    if(is.null(indicator.matrix)) return(NULL)
    
    if(nrow(indicator.matrix) > kpm.max.rows){
      showshinyalert(session, "kpm_status", paste("Error: Too many genes have been selected for this analysis. A maximum of ", kpm.max.rows, " is allowed with the current settings. Reduce the number of input genes by increasing the stringency of the previous analysis steps.", sep = ""), "danger")
      return(NULL)
    }
    
    list.of.ind.matrices <- list(indicator.matrix)
    
    #if(input$random.miRNA.test){
    #  showshinyalert(session, "kpm_status", "Generating matrices for random miRNAs", "info")    
    #  list.of.ind.matrices <- append(list.of.ind.matrices, list.of.random.mirna.indicator.matrices())
    #}      
    
    showshinyalert(session, "kpm_status", "Sending data to KeyPathwayMiner", "info")
    
    #if(input$kpm_ranged)
    #{
    #  Kmin <- input$kpm_lower_K
    #  Lmin <- input$kpm_lower_L
    #}
    #else{
      Kmin <- input$kpm_K
      
      if(input$screenType == "miRNA"){
        Lmin <- length(unique(mirna.targets()$mature_miRNA)) - input$kpm_L
      }
      else if(input$screenType == "siRNA")
      {
        Lmin <- 0
      }
      else if(input$screenType == "compound")
      {
        Lmin <- length(unique(drug.targets()$PubChem_CID)) - input$kpm_L
      }
      
    #}
    
    result <- foreach(ind.matrix = list.of.ind.matrices) %do%{
                    call.KPM(list(ind.matrix), 
                       url=keypathwayminer.url, 
                       ATTACHED_TO_ID=ATTACHED_TO_ID, 
                       strategy="INES",#input$kpm_strategy, 
                       algorithm="Greedy",#input$kpm_algorithm, 
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
                       computed.pathways=20)#input$kpm_pathways)
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
  paste("Click <a target='_blank' href='", kpm.url, ATTACHED_TO_ID, "&hideTitle=false", "'><u>here</u></a> to follow the progress of your run in KeyPathwayMiner web", sep="")
}

quest.result.url <- function(){
  kpm.url <- paste(keypathwayminer.url, "requests/quests?attachedToId=", sep="")
  paste("Click <a target='_blank' href='", kpm.url, ATTACHED_TO_ID, "&hideTitle=false", "'><u>here</u></a> to see the results of your run in KeyPathwayMiner web", sep="")  
}

check.KPM.result <- function(){

  #KPM hasn't been executed yet
  if(input$startKPMButton == 0){    
    return(FALSE)
  } 
  else if(is.null(currentQuest())) return(FALSE)
  
  #get current running status from KPM web
  quest.status <- getKpmRunStatus(keypathwayminer.url, currentQuest())
  if(is.null(quest.status)) return(FALSE)
  
  #update status
  showshinyalert(session, "kpm_status", paste("KeyPathwayMiner run is ", quest.status[["progress"]]*100, "% completed.", quest.progress.url()), "warning")  
  
  if(quest.status[["completed"]] == TRUE){
    #update quest.completed to avoid unneccessary polling    
    showshinyalert(session, "kpm_status", paste("KeyPathwayMiner run finished", quest.result.url()), "success")  
    
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


#modify hit list for KPM functionality if we are dealing with miRNAs
KPM.modify.hits <- reactive({
  hits <- outliers()
  
  if(input$accessionColType == "MI")
  {
    mimat <- as.data.frame(mirbaseMATURE)
    hits <- left_join(hits, mimat, by=c("mirna_id"))
  }
  return(hits)
})

#extract a single selected graph and process KPM data before plotting
kpm.graph.data <- reactive({
  hits <- KPM.modify.hits()
  kpm.res <- KPM.result()
  if(is.null(kpm.res)) return(NULL)
  
  selectedGraph <- NA
  if(!input$kpm_union_graph) selectedGraph <- input$kpm_selected_solution
  
  kpm.data <- NULL
  
  tryCatch(
    { 
      kpm.data <- prepare.kpm.output.for.plotting(kpm.res, KPM.indicator.matrix(), hits, input$screenType, selectedGraph)
    },
    error = function(e) { 
      showshinyalert(session, "kpm_status", e$message, "danger")
    }
  )
  return(kpm.data)
})

# Present nodes of currently selected graph as table
kpm.node.table <- reactive({
  graph.data <- kpm.graph.data()
  if(is.null(graph.data)) return(NULL)
  node.ids <- graph.data[[2]]
  node.ids <- node.ids[[1]]
  
  node.table <- dplyr::left_join(data.frame(name = node.ids, stringsAsFactors = FALSE), as.data.frame(org.Hs.egSYMBOL), by=c("name" = "gene_id"))
  
  if(input$screenType == "miRNA"){
    hit.list <- KPM.modify.hits() %>% dplyr::select(mature_name, category)
    node.table <- dplyr::left_join(node.table, hit.list, by=c("name"="mature_name"))
  }
  else if(input$screenType == "siRNA")
  {
    hit.list <- KPM.modify.hits() %>% dplyr::select(gene_id, category)
    node.table <- left_join(node.table, hit.list, by=c("name" = "gene_id"))
  }
  node.table[which(is.na(node.table$category)), "category"] <- "exception node"
  
  node.table <- arrange(node.table, category, name, symbol)
  return(node.table)
})