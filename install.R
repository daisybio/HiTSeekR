### installation ###
install.packages(c("shiny", "plyr", "dplyr", "ggplot2", 
                   "gplots", "scales", "gridExtra", 
                   "reshape2", "stringr", "VennDiagram", "qgraph",
                   "iterators", "foreach", "mirbase.db", "htmlwidgets",
                   "networkD3", "tidyr", "devtools", "XML", "R.utils"))
require(devtools)
install_github('ramnathv/rCharts')
install_github("AnalytixWare/ShinySky")

source("http://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db", ask=F)
biocLite("RmiR", ask=F)
biocLite("RTCA", ask=F)
biocLite("KEGG.db", ask=F)
biocLite("org.Dm.eg.db", ask=F)
biocLite("HTSanalyzeR", ask=F)
biocLite("KEGG.db", ask=F)
biocLite("org.Dm.eg.db", ask=F)
biocLite("HTSanalyzeR", ask=F)

### for parallelization ###
install.packages("doRedis")

#or
#install.packages("doParallel")
