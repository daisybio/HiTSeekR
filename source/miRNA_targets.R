load("data/microT.CDS.hsa.Rdata")

getTargets <- function(outliers, hits.min=1, at.least=2, 
                       group.miRNAs=T, group.miRNAs.threshold=2, 
                       get.gene.symbols=T, databases=NA)
{
  library(dplyr)
  outliers$mature_acc <- sub("mir", "miR", outliers$miRBase.ID.miRPlus.ID)
  
  targets <- semi_join(microT.hsa.processed, outliers, by="mature_acc")
  
  return(targets)
} 