require(shiny)
library(plyr)
library(dplyr)
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
require(shinysky)
library(iterators)
library(foreach)

source("source/heatmap.R")
source("source/RmiR2.R")
source("source/go.analysis.R")
source("source/Bscore.R")
source("source/posEffectNorm.R")
source("source/normalize.raw.data.R")
source("source/call_KPM.R")
source("source/highcharts_heatmap.R")
source("source/highcharts_scatterplot.R")
source("source/bayesian_hit_selection.R")
source("source/variance_based_hit_selection.R")
source("source/plot.miRNA.target.enrichment.graph.R")

options(shiny.maxRequestSize=30*1024^2)

# load mircancer database #
mircancer.database <- read.table("data/miRCancerSeptember2014.txt", sep="\t", header=T, quote="\"")

shinyServer(function(input, output, session) {
  
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
  
  source("ui/navbar_data_options.R", local = TRUE)
  source("ui/navbar_hits.R", local = TRUE)
  source("ui/navbar_consensus_hits.R", local = TRUE)
  source("ui/navbar_mirna_targets.R", local = TRUE)
  source("ui/gene_ontology.R", local = TRUE)
  source("ui/navbar_data.R", local = TRUE)
}) 
