#library(org.Hs.eg.db)
#library(RmiR)
#library(dplyr)
#library(foreach)

getTargets <- function(hits, rnah.pvalue.threshold=0.001, get.gene.symbols=F, databases=NA, diana.threshold=0.8){
  #repair ids
  
  hits <- na.omit(sub("mir", "miR", hits$mature_name))
  if("DIANA_microT_CDS" %in% databases){
    return(get_DIANA_microT_targets(hits, diana.threshold))
  }
  else if("DIANA_tarbase" %in% databases)
  {
    return(get_DIANA_tarbase6_targets(hits))
  }
  
  if("RNAhybrid_hsa" %in% databases) rnah.db <- src_sqlite(paste(data.folder, "rnahybrid.sqlite3", sep=""))
  else miRNA.db <- src_sqlite(RmiR.Hs.miRNA_dbfile())
  
  query.result <- foreach(db=databases, .combine=rbind) %do% {
    if(db == "RNAhybrid_hsa"){
      miRNA.targets <- tbl(rnah.db, "rnah")
      miRNA.targets <- filter(miRNA.targets, pvalue < rnah.pvalue.threshold)
    } 
    else miRNA.targets <- tbl(miRNA.db, db) 
    
    miRNA.targets <- miRNA.targets %>% dplyr::filter(mature_miRNA %in% hits) 
    miRNA.targets <- collect(miRNA.targets)
    miRNA.targets$gene_id <- as.integer(miRNA.targets$gene_id)
    
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

