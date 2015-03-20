#library(org.Hs.eg.db)
#library(RmiR)
#library(dplyr)
#library(foreach)

getTargets <- function(hits, rnah.pvalue.threshold=0.001, get.gene.symbols=F, databases=NA){
  #repair ids
  hits <- na.omit(sub("mir", "miR", hits$mature_name))
  
  if("RNAhybrid_hsa" %in% databases) rnah.db <- src_sqlite("data/rnahybrid.sqlite3")
  else miRNA.db <- src_sqlite(RmiR.Hs.miRNA_dbfile())
  
  query.result <- foreach(db=databases, .combine=rbind) %do% {
    if(db == "RNAhybrid_hsa"){
      miRNA.targets <- tbl(rnah.db, "rnah")
      miRNA.targets <- filter(miRNA.targets, pvalue < rnah.pvalue.threshold)
    } 
    else miRNA.targets <- tbl(miRNA.db, db) 
    
    miRNA.targets <- miRNA.targets %>% filter(mature_miRNA %in% hits) 
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
  
getRNAhybridTargets <- function(outliers, group.miRNAs=T, group.miRNAs.threshold=2, pvalue.threshold = 0.1){
    #repair ids
    outliers <- sub("mir", "miR", outliers$Sample)
    miRNA.db <- src_sqlite("rnah.sqlite3")
    miRNA.db <- tbl(miRNA.db, "rnah")
    miRNA.targets <- miRNA.db %>% filter(pvalue < pvalue.threshold)
    miRNA.targets <- miRNA.targets %>% filter(miRNA %in% outliers) 
    if(group.miRNAs){
            miRNA.targets <- miRNA.targets %>% 
                group_by(miRNA, Entrez) %>% 
                summarize(count=n(), source_db="RNAhybrid") %>% 
                filter(count > group.miRNAs.threshold)
        }
        miRNA.targets <- collect(miRNA.targets)
        
        return(miRNA.targets)
    
    
    #if(nrow(query.result > 20.000)) stop("Too many hits")
    #else 
    return(miRNA.targets)
}