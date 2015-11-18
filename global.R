demo.data.sets <- c("RNAi screen identifies Caspase4 as factor for TNFa signalling" = "TNFa_Casp4",
                    "miRNA mimics screen for vorinostat resistance genes" = "HCC_vorinostat_miRNA",
                    "RNAi screen for vorinostat resistance genes" = "HCC_vorinostat_siRNA",
                    "KRAS Synthetic Lethal Screen in HKE3" = "KRAS_synleth_compound")#,                    
                    #"BCSC/MaSC miRNA inhibitors" = "BCSC", 
                    #"Melanoma-Inhibiting miRNAs" = "A375_MTS", 
                    #"Drosophila M. Kc167 genome-wide siRNA screen"="DM_Kc167")

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
#library(org.Dm.eg.db)
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

### Additional shiny options ###
options(shiny.maxRequestSize=30*1024^2)

### Load mircancer database ###
mircancer.database <- read.table(paste(data.folder, "miRCancerSeptember2015.txt", sep=""), sep="\t", header=T, quote="\"")

### Load miRNA aliases. First try to read the most up to date version from mirbase website, use local copy as a fallback only ###
mirna.aliases <- NULL

tryCatch({
  con <- gzcon(url(paste("ftp://mirbase.org/pub/mirbase/CURRENT/aliases.txt.gz", sep="")))
  txt <- readLines(con)
  mirna.aliases <-  read.delim(textConnection(txt), header=F)
},
error = function(e) { mirna.aliases = read.delim(paste(data.folder, "aliases.txt", sep=""), header=F)
})