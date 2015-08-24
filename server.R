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
library(org.Dm.eg.db)
library(GO.db)
library(KEGG.db)
library(RCurl)

### Load all non-shiny source files ###
source("source/wellPositionToAlphaName.R")
source("source/heatmap.R")
source("source/RmiR.R")
source("source/go.analysis.R")
source("source/Bscore.R")
source("source/SSMD.R")
source("source/posEffectNorm.R")
source("source/normalize.raw.data.R")
source("source/call_KPM.R")
source("source/highcharts_heatmap.R")
source("source/highcharts_scatterplot.R")
source("source/bayesian_hit_selection.R")
source("source/variance_based_hit_selection.R")
source("source/plot.miRNA.target.enrichment.graph.R")
source("source/find.mimat.from.alias.R")
source("source/htsanalyzer.foreach.R")
source("source/htsanalyzer.gsea.optimized.R")
source("source/ggplot_smooth_func.R")

### Additional shiny options ###
options(shiny.maxRequestSize=30*1024^2)

### Load mircancer database #
mircancer.database <- read.table("data/mirCancerSeptember2014.txt", sep="\t", header=T, quote="\"")

### Shiny server ###
shinyServer(function(input, output, session) {
  
  ### Get parallel backend up and running ###
  
  if(require(doRedis)){  
    queueID <- paste(sample(c(LETTERS[1:6],0:9),8,replace=TRUE),collapse="")
    registerDoRedis(queueID)
    setChunkSize(value = 1000)
    startLocalWorkers(n=number.of.cores, queue=queueID)
  } else if(require(doParallel)){
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
  source("server/consensus_hits.R", local = TRUE)
    
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
  source("ui/navbar_consensus_hits.R", local = TRUE)
  source("ui/navbar_mirna_targets.R", local = TRUE)
  source("ui/navbar_drug_targets.R", local = TRUE)
  source("ui/navbar_data.R", local = TRUE)
  source("ui/navbar_gene_set.R", local = TRUE)
  source("ui/navbar_enable_tabs.R", local = TRUE)
  
  shinyjs::disable("qc")
  
  ### Cleanup, close parallel backend clusters if necessary ###
  cancel.onSessionEnded <- session$onSessionEnded(function() {    
    if(require(doRedis)) removeQueue(queueID)
    else if(require(doParallel)) stopCluster(cl)
  })
}) 
