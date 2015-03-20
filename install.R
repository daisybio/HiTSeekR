### installation ###

install.packages(c("shiny", "plyr", "dplyr", "ggplot2", 
                   "gplots", "scales", "gridExtra", 
                   "reshape2", "stringr", "VennDiagram", "qgraph",
                   "iterators", "foreach", "mirbase.db", "htmlwidgets",
                   "networkD3", "tidyr", "devtools"))
require(devtools)
install_github('rCharts', 'ramnathv')
install_github("AnalytixWare/ShinySky")

source("http://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
biocLite("RmiR")
biocLite("RTCA")

### for parallelization ###
install.packages("doRedis")

#or
#install.packages("doParallel")
