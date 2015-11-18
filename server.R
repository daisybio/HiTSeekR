### Shiny server ###
shinyServer(function(input, output, session) {
  
  ### Start page ###
  source("ui/frontpage.R", local = TRUE)
  
  screenType <- reactiveValues(type = NULL)  
  
  observeEvent(input$siRNA, {
    updateSelectInput(session, "screenType", "Type of screen", c("Gene silencing" = "siRNA", "miRNA inhibitor / mimics" = "miRNA", "Compound screen" = "compound"), "siRNA")
    updateSelectInput(session, "dataset", "Select a demo dataset", choices = c("none selected" = "none selected", demo.data.sets[c(1,3)]), "none selected")
    session$sendCustomMessage(type = "disableScreenType", message=list(empty=""))
  })
  
  observeEvent(input$miRNA, {
    updateSelectInput(session, "screenType", "Type of screen", c("Gene silencing" = "siRNA", "miRNA inhibitor / mimics" = "miRNA", "Compound screen" = "compound"), "miRNA")
    updateSelectInput(session, "dataset", "Select a demo dataset", choices = c("none selected" = "none selected", demo.data.sets[2]), "none selected")    
    session$sendCustomMessage(type = "disableScreenType", message=list(empty=""))
  })  
  
  observeEvent(input$compound, {
    updateSelectInput(session, "screenType", "Type of screen", c("Gene silencing" = "siRNA", "miRNA inhibitor / mimics" = "miRNA", "Compound screen" = "compound"), "compound")
    updateSelectInput(session, "dataset", "Select a demo dataset", choices = c("none selected" = "none selected", demo.data.sets[4]), "none selected")
    session$sendCustomMessage(type = "disableScreenType", message=list(empty=""))
  })  
  
  output$uiOutput_frontpage <- renderUI({ 
    #prepare progress bar
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Loading HiTSeekR...")    
    if(!is.null(screenType)) do.call(navbarPage, elts)
    else return(NULL)
  })  
  
  ### Get parallel backend up and running ###
  
  if(use.redis && require(doRedis)){  
    queueID <- paste(sample(c(LETTERS[1:6],0:9),8,replace=TRUE),collapse="")
    registerDoRedis(queueID, host=redis.host, nodelay=FALSE)
    setChunkSize(value = 50)
    startLocalWorkers(n=number.of.cores, queue=queueID, host=redis.host, timeout=2, nodelay=FALSE)
  } else if(require(doParallel)){
    if(number.of.cores == "auto") number.of.cores <- parallel::detectCores() 
    if(number.of.cores < 1) number.of.cores = 1
    cl <- makeCluster(number.of.cores)
    registerDoParallel(cl)
  }
    
  ### DATA PROCESSING ###
  
  #load the data
  source("server/load_data.R", local = TRUE)  
  
  #process data
  source("server/process_data.R", local=TRUE)
  
  #data columns
  dataColumns <- reactive({
    colnames(rawData())
  })
  
  #default data options for demo data sets
  source("server/default_options.R", local = TRUE)
    
  # filter and summarize
  source("server/filter_data.R", local = TRUE)
  
  #find hits
  source("server/find_hits.R", local = TRUE)
  
  #consensus hits
  #source("server/consensus_hits.R", local = TRUE)
    
  #miRNA targets
  source("server/mirna_targets.R", local = TRUE)
  
  #drug targets
  source("server/drug_targets.R", local = TRUE)
  
  #KeyPathwayMiner
  source("server/KPM.R", local = TRUE)
  
  ## mirCancerDB ##
  source("server/mircancer.R", local = TRUE)
  
  ## HTSanalyzer ##
  source("server/htsanalyzer.R", local = TRUE)
  
  ### Plots ###
  
  source("server/output/mRNA_miRNA_interaction_graph.R", local = TRUE)  
  source("server/output/plots.R", local = TRUE)
  
  #### Tables ###
  
  source("server/output/tables.R", local = TRUE)
  
  ### Downloads ###

  source("server/output/downloads.R", local = TRUE)

  ### User Interface ###
  
  source("ui/navbar_data_options.R", local = TRUE)
  source("ui/navbar_quality_control.R", local = TRUE)
  source("ui/navbar_hits.R", local = TRUE)
  #source("ui/navbar_consensus_hits.R", local = TRUE)
  source("ui/navbar_mirna_targets.R", local = TRUE)
  source("ui/navbar_drug_targets.R", local = TRUE)
  source("ui/navbar_data.R", local = TRUE)
  source("ui/navbar_gene_set.R", local = TRUE)
  source("ui/navbar_enable_tabs.R", local = TRUE)
  
  ### Load help pages ###
  
  source("server/output/help.R", local = TRUE)
  
  ### Cleanup, close parallel backend clusters if necessary ###
  cancel.onSessionEnded <- session$onSessionEnded(function() {    
    if(use.redis && require(doRedis)) removeQueue(queueID)
    else if(require(doParallel)) stopCluster(cl)
  })
}) 
