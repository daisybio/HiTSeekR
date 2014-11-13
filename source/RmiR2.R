library(org.Hs.eg.db)
library(RmiR)
library(dplyr)
library(foreach)

getTargets <- function(outliers, group.miRNAs=T, group.miRNAs.threshold=2, get.gene.symbols=T, databases=NA){
  #repair ids
  outliers <- sub("mir", "miR", outliers$miRBase.ID.miRPlus.ID)
  miRNA.db <- src_sqlite(RmiR.Hs.miRNA_dbfile())
  
  query.result <- foreach(db=databases, .combine=rbind) %do% {
    miRNA.targets <- tbl(miRNA.db, db) %>% filter(mature_miRNA %in% outliers) 
    if(group.miRNAs){
      miRNA.targets <- miRNA.targets %>% 
      group_by(mature_miRNA, gene_id) %>% 
      summarize(count=n(), source_db=db) %>% 
      filter(count > group.miRNAs.threshold)
    }
    miRNA.targets <- collect(miRNA.targets)
    #get gene symbols and add them to the result table
    if(nrow(miRNA.targets) > 0 & get.gene.symbols){      
      gene_ids <- unique(miRNA.targets$gene_id)
      gene_symbols <- as.character(mget(as.character(gene_ids), org.Hs.egSYMBOL, ifnotfound=NA))                                   
      names(gene_symbols) <- gene_ids
      miRNA.targets$gene_symbol <- gene_symbols[miRNA.targets$gene_id]
    }
    
    return(miRNA.targets)
  }
  
  #if(nrow(query.result > 20.000)) stop("Too many hits")
  #else 
  return(query.result)
}
  