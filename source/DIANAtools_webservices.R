get_DIANA_tarbase6_targets <- function(miRNAs)
{
  #read for each microRNA targets as XML from DIANA web service
  result <- foreach(miRNA = miRNAs, .combine=rbind) %do%
  {
    xmlDOC <- xmlTreeParse(getURL(paste("http://62.217.127.8/DianaTools/tarbaseApi?mirnas=", miRNA, sep="")))
    xmlDOM <- xmlRoot(xmlDOC)
    
    #check if there are results
    if(as.integer(xmlValue(xmlDOM[["info"]][["result-size"]])) == 0) return(NULL)
    
    #extract target genes from XML
    result <- foreach(target = 1:length(xmlDOM[["results"]]), .combine=rbind) %do%
    {
      xmlAttrs(xmlDOM[["results"]][[target]])
    }
    return(result)
  }
  
  #add entrez id
  result <- as.data.frame(result)
  row.names(result) <- NULL
  result$Ensembl <- NA
  result[which(grepl("ENSG[0-9]", result$geneId)), "Ensembl"] <- grep("ENSG[0-9]", result$geneId, value=TRUE)
  result <- dplyr::left_join(result, as.data.frame(org.Hs.egENSEMBL), by=c("Ensembl" = "ensembl_id")) %>%
    dplyr::rename(gene_symbol=geneName, gene_name_in_publication=geneId, mature_miRNA = mirnaName) 
  
  return(result)
}

get_DIANA_microT_targets <- function(miRNAs, threshold = 0.8)
{
  #read for each microRNA targets as XML from DIANA web service
  result <- foreach(miRNA = miRNAs, .combine=rbind) %do%
  {
    xmlDOC <- xmlTreeParse(getURL(paste("http://62.217.127.8/DianaTools/microT_CDSApi?mirnas=", miRNA, "&threshold=", threshold, sep="")))
    xmlDOM <- xmlRoot(xmlDOC)
    
    #extract target genes from XML
    result <- foreach(target = 1:length(xmlDOM[["results"]]), .combine=rbind) %do%
    {
      xmlAttrs(xmlDOM[["results"]][[target]])
    }
    return(result)
  }
  
  #split ensembl ID and gene symbol
  result <- as.data.frame(result) %>% separate(geneName, c("Ensembl", "GeneSymbol"), " ") %>% dplyr::mutate(GeneSymbol = str_sub(GeneSymbol, 2, -2))
  
  #add entrez id
  result <- dplyr::left_join(result, as.data.frame(org.Hs.egENSEMBL), by=c("Ensembl" = "ensembl_id")) %>%
    dplyr::select(mature_miRNA = mirnaName, Ensembl, GeneSymbol, gene_id, score) %>% dplyr::arrange(mature_miRNA, score)
  
  return(result)
}

get_DIANA_mirPath <- function(miRNAs, threshold = 0.8, geneIntersectionCutoff=2, selection=1, fdr=0.05, conservative=TRUE)
{
  if(is.null(miRNAs)) return(NULL)
  miRNAs <- paste(na.omit(miRNAs$mature_name), collapse=" ")
  miRNAs.encoded <- str_replace_all(miRNAs, " ", "%20")
  xmlDOC <- xmlTreeParse(getURL(paste("http://62.217.127.8/DianaTools/mirpathApi?mirnas=", miRNAs.encoded, 
                                      "&methods=microT-CDS&stats=", conservative, 
                                      "&species=human",
                                      "&false_rate=", fdr, "&selection=", selection, "&inter_mir=", geneIntersectionCutoff, 
                                      "&threshold=", threshold, sep="")))
  xmlDOM <- xmlRoot(xmlDOC)
  
  #extract kegg pathways from XML
  result <- foreach(record = 1:length(xmlDOM[["results"]]), .combine=rbind) %do%
  {
    entry <- xmlDOM[["results"]][[record]]
    entry.miRNAs <- foreach(miR = 1:length(entry), .combine=paste) %do% {
      sub("\\|.*", "", xmlAttrs(entry[[miR]])[[1]])
    }
    entry.miRNAs <- str_replace_all(entry.miRNAs, " ", "\n")
    entry <- xmlAttrs(entry)
    entry <- c(entry, miRNAs = entry.miRNAs)
    return(entry)
  }

  return(result)
}
