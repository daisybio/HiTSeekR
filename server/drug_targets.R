drug.targets <- reactive({
  
  exp.data <- data()
  stitch.db <- src_sqlite("data/stitch_hsa_protein_chemical_links_v4.0.sqlite3")
  targets <- tbl(stitch.db, "hs")
  
  targets %>% filter(PubChemID %in% exp.data$Accession)
})