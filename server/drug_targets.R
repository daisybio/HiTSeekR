selectedHitList <- reactive({
  if(input$drugUseConsensus == "hit list") data <- outliers()
  else data <- consensusHitList()
  return(data)
})

convertToCid <- function(data, type){
  if(type == "PubChem_CID"){
    data$PubChem_CID <- data$Accession
    return(data)
  }
  
  compound.mapping.db <- src_sqlite("data//compound.mapping.sqlite3")
  
  if(type == "Chembank")
  {
      compound.mapping <- tbl(compound.mapping.db, "chembank_to_cid")
      result <- filter(compound.mapping, Chembank_ID %in% data$Accession)
      result <- as.data.frame(result) %>% filter(complete.cases(.))
      result <- left_join(data, result, by=c("Accession" = "ChemBank_ID"))
      return(result)
  }
  else if(type == "PubChem_SID")
  {
      compound.mapping <- tbl(compound.mapping.db, "chembank_to_sid_cid")
      result <- filter(compound.mapping, PubChem_SID %in% data$Accession)
      result <- as.data.frame(result) %>% filter(complete.cases(.)) %>% select(-ChemBank_ID)
      result <- left_join(data, result, by=c("Accession" = "PubChem_SID"))    
      return(result)
  }
  
  stop("Do not know how to handle unknown compound / drug ID type.")
}

drug.targets <- reactive({
  
  hits <- selectedHitList()
  
  #convert ids to pubchem compound ids (CID) as used in STITCH
  hits <- convertToCid(hits, input$accessionColType)
  
  #add CID prefix and leading zeros for STITCH db
  hits[!is.na(hits$PubChem_CID), "PubChem_CID"] <- paste("CID", formatC(as.integer(hits[!is.na(hits$PubChem_CID), "PubChem_CID"]), width=9, flag="0"), sep="")
  
  tryCatch({
  stitch.db <- src_sqlite("data/stitch_hsa_protein_chemical_links_v4.0.sqlite3")
  targets <- tbl(stitch.db, "hs")
  hits.wo.nas <- na.omit(hits$PubChem_CID)
  
  targets <- filter(targets, PubChemID %in% hits.wo.nas)
  hits <- hits %>% dplyr::select(Sample, PubChem_CID) 
  hits <- left_join(hits, targets, by=c("PubChem_CID" = "PubChemID"), copy=T)
  }, error = function(e){ showshinyalert(session, "general_status", paste("error:", e$message), "danger") })
  
  if(nrow(hits) == 0) stop("No target has been found. Try to reduce stringency to increase number of hits.")
  
  return(hits)
})
