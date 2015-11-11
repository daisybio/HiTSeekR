getRNAhybridTargetCounts <- function(genes=NULL, pval.threshold)
{
  rnah.db <- src_sqlite(paste(data.folder, "rnahybrid.sqlite3", sep=""))
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