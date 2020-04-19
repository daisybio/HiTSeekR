### installation ###

# this file is thought as a guide for users to install all necessary packages to run HiTSeekR locally.

#CRAN packages
install.packages(c("shiny", "plyr", "dplyr", "ggplot2", 
                   "gplots", "scales", "gridExtra", "lazyeval",
                   "reshape2", "stringr", "VennDiagram", "qgraph",
                   "iterators", "foreach", "htmlwidgets", "DBI",
                   "networkD3", "tidyr", "devtools", "XML", "R.utils", "BiocManager"))

#bioconductor packages
BiocManager::install(c("org.Hs.eg.db", "RmiR", "RTCA", "KEGG.db", "mirbase.db", "reactome.db", "GO.db"))

#github packages
library(devtools)
install_github('ramnathv/rCharts')
install_github("AnalytixWare/ShinySky")

#for gene set analysis
#biocLite("HTSanalyzeR", ask=F)

#load github version that has been modified for foreach support and shiny progress bars
install_github("NanoCAN/HTSanalyzeR")

### for parallelization, recommended to use one of the two solutions ###

#NOTE: using redis as parallel backend means a redis service has to be running and should be configured in config.R
#The advantage of using redis is that the progress bars are dynamically updated for long running computations.

#install.packages("doRedis")
install_github("bwlewis/doRedis")

#or as a fallback
#install.packages("doParallel")
