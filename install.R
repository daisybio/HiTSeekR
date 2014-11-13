#these are the required packages for running RNAice
install.packages(c("shiny", "qgraph", "ggplot2","grid","gridExtra", "scales", "devtools", "data.table", "DBI", "gplots", "VennDiagram"))

library(devtools)
install_github('rCharts', 'ramnathv', ref='dev')

source("http://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
biocLite("RmiR")
