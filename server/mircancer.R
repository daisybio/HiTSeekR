# find entries for currently selected hits #
hits.mircancer <- reactive({
  hits <- outliers()
  hits <- hits[,c("mirna_id", "category", "Experiment", "Accession", "Sample")]
  mirdb <- mircancer.database
  showshinyalert(session, "mircancer_status", paste("This information was extracted from the",
    mircancer.version, "release of <a href='http://mircancer.ecu.edu/'>miRCancer</a>. To find the corresponding PubMed articles browse the miRNA in question by clicking on its id."),"info")
  return(as.data.frame(left_join(hits, mirdb, by = c("mirna_id" = "mirId"))))
})