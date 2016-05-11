### Configuration ###
source("config.R")

### Load required packages ###
library(shiny)
library(plyr)
library(dplyr)
library(rCharts)
library(gplots)
library(ggplot2)
library(grid)
library(gridExtra)
library(scales)
library(reshape2)
library(RmiR)
library(stringr)
library(VennDiagram)
library(qgraph)
library(RTCA)
library(shinysky)
library(iterators)
library(foreach)
library(mirbase.db)
library(htmlwidgets)
library(networkD3)
library(tidyr)
library(org.Hs.eg.db)
library(HTSanalyzeR)
library(GSEABase)
library(GO.db)
library(KEGG.db)
library(RCurl)
library(lazyeval)
library(XML)

### Load all non-shiny source files ###
source("source/wellPositionToAlphaName.R")
source("source/heatmap.R")
source("source/go.analysis.R")
source("source/Bscore.R")
source("source/SSMD.R")
source("source/posEffectNorm.R")
source("source/find_hits_call.R")
source("source/normalize.raw.data.R")
source("source/call_KPM.R")
source("source/highcharts_heatmap.R")
source("source/highcharts_scatterplot.R")
source("source/bayesian_hit_selection.R")
source("source/variance_based_hit_selection.R")
source("source/plot.kpm.result.graph.R")
source("source/find.mimat.from.alias.R")
source("source/htsanalyzer.reactomeGeneSets.R")
source("source/ggplot_smooth_func.R")
source("source/DIANAtools_webservices.R")
source("source/RNAhybrid.R")
source("source/RmiR.R")
devtools::source_url("https://raw.githubusercontent.com/NanoCAN/SAVANAH/master/R/import.R")

### Additional shiny options ###
options(shiny.maxRequestSize=30*1024^2)

### Load mircancer database ###
mircancer.database <- read.table(paste(data.folder, "miRCancerSeptember2015.txt", sep=""), sep="\t", header=T, quote="\"")

### Load miRNA aliases. First try to read the most up to date version from mirbase website, use local copy as a fallback only ###
mirna.aliases <- NULL

### Load miRNA gene target significance counts for RmiR ###
load("data/rmir.counts.RData")

tryCatch({
  con <- gzcon(url(paste("ftp://mirbase.org/pub/mirbase/CURRENT/aliases.txt.gz", sep="")))
  txt <- readLines(con)
  mirna.aliases <-  read.delim(textConnection(txt), header=F)
},
error = function(e) { mirna.aliases = read.delim(paste(data.folder, "aliases.txt", sep=""), header=F)
})

### Shiny server ###
shinyServer(function(input, output, session) {
  
  ### Check if we are called by SAVANAH ###
  # Parse the GET query string
  output$queryText <- renderText({
    query <- parseQueryString(session$clientData$url_search)
    
    # Return a string with key-value pairs
    paste(names(query), query, sep = "=", collapse=", ")
  })
  
  isRemoteCall <- reactive({
    query <- parseQueryString(session$clientData$url_search)
    if(is.null(query$plateSecurityTokens) || is.null(query$screenType)) return(FALSE)
    else return(TRUE)
  })
  
  getBaseUrl <- reactive({
    query <- parseQueryString(session$clientData$url_search)
    if(length(query$baseUrl > 0)) baseUrl <- query$baseUrl
    else baseUrl <- "http://localhost:8080/SAVANAH/readoutExport/"
    return(baseUrl)
  })
  
  loadAllReadouts <- reactive({ 
    baseUrl <- getBaseUrl()

    query <- parseQueryString(session$clientData$url_search)    
    if(query$plateSecurityTokens != ""){
      plateTokens <- str_split(query$plateSecurityTokens, "\\|")[[1]]
    } 
    else return(NULL)
    
    batch.import.readouts(plateSecurityTokens=plateTokens,baseUrl = baseUrl)
  })
  
  
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
  
  observeEvent(input$screenType, {
    session$sendCustomMessage(type = "disableLoadingScreen", message=list(empty=""))
    
    if(isRemoteCall()){
      query <- parseQueryString(session$clientData$url_search)
      updateSelectInput(session, "screenType", "Type of screen", c("Gene silencing" = "siRNA", "miRNA inhibitor / mimics" = "miRNA", "Compound screen" = "compound"), query$screenType)
      updateSelectInput(session, "dataset", "Select a demo dataset", choices = c("none selected" = "none selected", demo.data.sets[c(1,3)]), "none selected")
      session$sendCustomMessage(type = "disableScreenType", message=list(empty=""))
      #session$sendCustomMessage(type = "disableInputTab", message=list(empty=""))
      session$sendCustomMessage(type = "disableNavTab", "1")
      session$sendCustomMessage(type = "enableNavTab", "2")
      session$sendCustomMessage(type = "enableNavTab", "3")
      session$sendCustomMessage(type = "enableNavTab", "4")
      #get processedData filled up correctly
    }
  })
  
  output$uiOutput_frontpage <- renderUI({ 
    #prepare progress bar
    progress <- shiny::Progress$new()
    progress$set(message = "Loading HiTSeekR...")   
    
    if(!is.null(screenType)){
      do.call(navbarPage, elts)
    }
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
