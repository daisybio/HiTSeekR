### installation ###
install.packages(c("shiny", "plyr", "dplyr", "ggplot2", 
                   "gplots", "scales", "gridExtra", 
                   "reshape2", "stringr", "VennDiagram", "qgraph",
                   "iterators", "foreach", "htmlwidgets", "shinyjs",
                   "networkD3", "tidyr", "devtools", "XML", "R.utils"))

library(devtools)
install_github('ramnathv/rCharts')
install_github("AnalytixWare/ShinySky")

source("http://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db", ask=F)
biocLite("RmiR", ask=F)
biocLite("RTCA", ask=F)
biocLite("KEGG.db", ask=F)
biocLite("org.Dm.eg.db", ask=F)
biocLite("mirbase.db")

#for gene set analysis
#biocLite("HTSanalyzeR", ask=F)

#load github version that has been modified for foreach support and shiny progress bars
install_github("NanoCAN/HTSanalyzeR")

### for parallelization ###
#install.packages("doRedis")
install_github("bwlewis/doRedis")

#or
#install.packages("doParallel")
