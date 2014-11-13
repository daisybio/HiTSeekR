require(shiny)
require(plyr)
require(rCharts)
require(gplots)
require(ggplot2)
require(grid)
require(gridExtra)
require(scales)
require(reshape2)
require(RmiR)
require(stringr)
require(VennDiagram)
require(qgraph)
require(RTCA)

source("source/heatmap.R")
source("source/RmiR2.R")
source("source/go.analysis.R")
source("source/Bscore.R")
source("source/posEffectNorm.R")
source("source/normalize.raw.data.R")
source("source/call_KPM.R")

options(shiny.maxRequestSize=30*1024^2)

shinyServer(function(input, output) {
  
  ### DATA PROCESSING ###
  
  #load the data
  source("server/load_data.R", local = TRUE)  
  
  #data columns
  dataColumns <- reactive({
    colnames(rawData())
  })
  
  #data option default
  dataOptionDefaults <- reactive({
    if(input$dataset == "BCSC"){
      return(c("screenType" = "miRNA",
        "sampleCol" = "Sample", 
        "posColType" = "alpha",
        "posCol" = "Well.position",
        "accColType" = "MIMAT",
        "accCol" = "miRBase.accession",
        "measurementCol" = "CTB",
        "replicateCol" = "Replicate",
        "plateCol" = "Plate",
        "expCol" = "cellType",
        "ctrlCol" = "Control"
      ))
    }
    else return(NULL)
  })
  
  #process data
  processedData <- reactive({
    data <- rawData()
    sampleCol <- data[,input$sampleCol]
    if(input$positionColType == "alpha")
    {
      wellAlpha <- repairAlphaName(data[,input$positionCol])
      rowCol <- alphaNames2Pos(wellAlpha)
    }
    else if(input$positionColType == "numeric")
    {
      #TODO
    }
    replicateCol <- data[,input$replicateCol]
    measurementCol <- data[,input$measurementCol]
    accessionCol <- data[,input$accessionCol]
    plateCol <- data[,input$plateCol]
    experimentCol <- data[,input$experimentCol]
    controlCol <- data[,input$controlCol]
    
    processedData <- cbind(experimentCol, sampleCol, accessionCol, plateCol, wellAlpha, rowCol, replicateCol, controlCol, measurementCol)
    colnames(processedData) <- c("Experiment", "Sample", "Accession", "Plate", "Well.position", "Row", "Column", "Replicate", "Control", "Raw")
    
    if(input$hasControls && !(is.null(input$controlCol)))
    {
      normalizeRawData(processedData, control.based=T, pos.ctrl=input$posCtrl, neg.ctrl=input$negCtrl)
    }
    else{
      normalizeRawData(processedData, control.based=F)    
    }    
  })
  
  #identifier column types
  accTypes <- reactive({
    #if(input$screenType == "siRNA") return(c("RefSeq", "Entrez", "GeneSymbol"))
    #else if(input$screenType == "miRNA") return(c("MI", "MIMAT"))
    return(c("RefSeq", "Entrez", "MIMAT", "MI"))
  })
  
  output$uiOutput_data_options <- renderUI({
    elements <- list(
      selectInput("screenType", "Type of screen", c("Gene knockout (e.g. siRNA)" = "siRNA", "miRNA inhibitor / mimics" = "miRNA"), dataOptionDefaults()[["screenType"]]),
      selectInput("sampleCol", "Sample Name Column", dataColumns(), dataOptionDefaults()[["sampleCol"]]),
      selectInput("plateCol", "Plate Column", dataColumns(), dataOptionDefaults()[["plateCol"]]),
      selectInput("positionColType", "Position Column Type", c("Alpha well names" = "alpha", "Numeric" = "numeric"), dataOptionDefaults()[["posColType"]]),
      selectInput("positionCol", "Position Column", dataColumns(), dataOptionDefaults()[["posCol"]]),
      selectInput("accessionColType", "Accession Column Type", accTypes(), dataOptionDefaults()[["accColType"]]),
      selectInput("accessionCol", "Accession Column", dataColumns(), dataOptionDefaults()[["accCol"]]),
      selectInput("measurementCol", "Measurement Column", dataColumns(), dataOptionDefaults()[["measurementCol"]]),
      selectInput("replicateCol", "Replicate Column", dataColumns(), dataOptionDefaults()[["replicateCol"]]),
      selectInput("experimentCol", "Experiment Column", dataColumns(), dataOptionDefaults()[["expCol"]]),
      checkboxInput("hasControls", "Are controls included?", TRUE),
      conditionalPanel(condition = "input.hasControls",
        selectInput("controlCol", "Control Column", dataColumns(), dataOptionDefaults()[["ctrlCol"]])),        
        uiOutput("uiOutput_controls")
    )
    do.call(wellPanel, elements)
  })
  
  output$uiOutput_controls <- renderUI({
    if(is.null(input$controlCol) || !input$hasControls) return(NULL)
    
    data <- rawData()
    
    ctrlTypes <- unique(as.character(data[,input$controlCol]))
    
    ctrlSelects <- list(
    selectInput("posCtrl", "Positive Control", ctrlTypes, multiple=T),
    selectInput("negCtrl", "Negative Control", ctrlTypes, multiple=T)
    )
    do.call(wellPanel, ctrlSelects)
  })
  
  # filter and summarize
  source("server/filter_data.R", local = TRUE)
  
  #find hits
  source("server/find_hits.R", local = TRUE)
  
  #consensus hits
  source("server/consensus_hits.R", local = TRUE)
    
  #miRNA targets
  source("server/mirna_targets.R", local = TRUE)
  
  #KeyPathwayMiner
  source("server/KPM.R", local = TRUE)
  
  ## mirCancerDB ##
  source("server/mircancer.R", local = TRUE)
  
  ## gene ontology enrichment analysis ##
  goEnrichment <- reactive({
    goEnrichmentAnalysis(targets(), goUpOrDown=input$goUpOrDown, goDomain=input$goDomain, goScoringThreshold=input$goSignThreshold, orderMethod=input$goOrderMethod, topNodes=input$goTopNodes)
  })
  
  ### Plots ###
  
  source("server/output/mRNA_miRNA_interaction_graph.R", local = TRUE)  
  source("server/output/plots.R", local = TRUE)
  
  #### Tables ###
  
  source("server/output/tables.R", local = TRUE)
  
  ### Downloads ###

  source("server/output/downloads.R", local = TRUE)

  ### User Interface ###
  
  source("ui/navbar_hits.R", local = TRUE)
  source("ui/navbar_consensus_hits.R", local = TRUE)
  source("ui/navbar_mirna_targets.R", local = TRUE)
  source("ui/gene_ontology.R", local = TRUE)
  source("ui/navbar_mircancer.R", local = TRUE)
  source("ui/navbar_data.R", local = TRUE)
}) 
