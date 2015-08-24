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

getRNAhybridTargetCounts <- function(genes=NULL, pval.threshold)
{
  rnah.db <- src_sqlite("data/rnahybrid.sqlite3")
  dbname <- switch(as.character(pval.threshold), 
                   "0.05" = "mircounts05",
                   "0.01" = "mircounts01",
                   "0.001" = "mircounts001",
                   "1e-04" = "mircounts0001")
  mirna.counts <- tbl(rnah.db, dbname)
  if(is.null(genes)) return(dplyr::select(mirna.counts, gene))
  else return(mirna.counts %>% filter(gene %in% genes))  
}

getRNAhybridNumOfmiRNAs <- function(){
assign('n_distinct', function(x) {build_sql("COUNT(DISTINCT ", x, ")")}, envir=base_agg)
rnah.db <- src_sqlite("data/rnahybrid.sqlite3")
mirna.targets <- tbl(rnah.db, "rnah")
result <- as.data.frame(mirna.targets %>% dplyr::summarise(n_distinct(mature_miRNA)))
return(result[1,1])
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