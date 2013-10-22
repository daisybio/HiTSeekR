#function to merge column name into one
my.aggregate <- function(x)
{
  paste(unique(setdiff(x,NA)), collapse="/")
}

getTargets <- function(outliers, hits.min=1, at.least=2, group.miRNAs=T, group.miRNAs.threshold=2, get.gene.symbols=T, databases=NA){
  require(RmiR)
  require(plyr)
  miRNA.targets <- data.frame()
  #repair ids
  outliers$miRBase.ID.miRPlus.ID <- sub("mir", "miR", outliers$miRBase.ID.miRPlus.ID)
  
  if(length(databases) == 1) databases <- dbListTables(RmiR.Hs.miRNA_dbconn())
  
  #query all dbs
  mirnas <- unique(outliers$miRBase.ID.miRPlus.ID)
  firstLoop <- TRUE

  #group miRNAs together
  if(group.miRNAs) query <- "SELECT gene_id, COUNT(DISTINCT mature_miRNA) as miRNA_count, REPLACE(GROUP_CONCAT(DISTINCT mature_miRNA), ',', '/') AS miRNA_list, SUM(sum_hits) AS total_hits, MIN(db_count) AS min_db_count FROM ("
  
  else query <- ""
  
  #group miRNA / gene pairs and sum up number of hits in various databases, as well as number of databases where this interaction is reported
  query <- paste(query, "SELECT mature_miRNA, gene_id, SUM(hits) AS sum_hits, COUNT(DISTINCT database) AS db_count, REPLACE(GROUP_CONCAT(DISTINCT database), ',', '/') AS db_list FROM (", sep="")
  
  #query each database for miRNA targets, merge results using UNION and by adding a database column on the fly
  for(database in databases)
  {
    if(!firstLoop){ query <- paste(query, " UNION ALL ", sep="") }
    firstLoop <- FALSE
    query <- paste(query, "SELECT mature_miRNA, CAST(gene_id AS NUMBER) AS gene_id, COUNT(gene_id) AS hits, '", database, 
                   "' AS database FROM ", database, " WHERE mature_miRNA IN ", 
                   paste("('", paste(mirnas, collapse="','"), "')", sep="")  ,
                   " GROUP BY mature_miRNA, gene_id", sep="")
  }
  
  #finish grouping SQL statement for gene/miRNA pairs, filter for minimum number of databases and hits
  query <- paste(query, ") GROUP BY gene_id, mature_miRNA HAVING db_count >= ", at.least, " AND sum_hits >= ", hits.min, sep="") 
  
  #finish grouping SQL statement for grouping miRNAs together for each gene, filter for minimum number of interactions found per gene
  if(group.miRNAs) query <- paste(query, ") GROUP BY gene_id HAVING COUNT(DISTINCT mature_miRNA) >= ", group.miRNAs.threshold,sep="")
  
  #execute the final SQL statement
  miRNA.targets <- dbGetQuery(RmiR.Hs.miRNA_dbconn(), query)
  
  #get gene symbols and add them to the result table
  if(nrow(miRNA.targets) > 0 && get.gene.symbols){
    require(org.Hs.eg.db)
    miRNA.targets$gene_symbol <- paste("<a href='http://www.ncbi.nlm.nih.gov/gene/?term=", 
                                       miRNA.targets$gene_id, "%5Buid%5D", "' target='_blank'>",
      as.character(mget(as.character(miRNA.targets$gene_id), org.Hs.egSYMBOL, ifnotfound=NA)),
                                       "</a>", sep="")
  }
  #return result
  return(miRNA.targets)
}