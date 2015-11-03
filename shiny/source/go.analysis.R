#go.targets <- read.csv("targets_miranda_mirbase_mirtarget2_pictar_tarbase_targetscan_1.5_SD_ (6).csv", header=T)

goEnrichmentAnalysis <- function(go.targets, goUpOrDown="both", goDomain="BP", goScoringThreshold=3, orderMethod="elim.KS", topNodes=500){
  require(topGO)
  require(Rgraphviz)
  
  #filter genes (promotors, suppressors or both)
  if(goUpOrDown == "promotors")
  {
    geneList <- abs(go.targets[go.targets$promotor_vs_suppressor > 0, "promotor_vs_suppressor"])
    names(geneList) <- go.targets[go.targets$promotor_vs_suppressor > 0, "gene_id"]
  } else if(goUpOrDown == "suppressors")
  {
    geneList <- abs(go.targets[go.targets$promotor_vs_suppressor < 0, "promotor_vs_suppressor"])
    names(geneList) <- go.targets[go.targets$promotor_vs_suppressor < 0, "gene_id"]
  } else if(goUpOrDown == "both")
  {
    geneList <- abs(go.targets$promotor_vs_suppressor)
    names(geneList) <- go.targets$gene_id
  }
  
  if(length(geneList) == 0) return(data.frame(error="Target list does not contain any genes for this selection"))
  
  #filter function for significant genes
  topGenes <- function(geneList)
  {
    return(geneList > goScoringThreshold)
  }
  
  #create topGO data object
  go.topGO <- new("topGOdata", ontology=goDomain, allGenes=geneList, geneSel = topGenes, annot=annFUN.org, mapping="org.Hs.eg.db")
  
  #perform statistical tests
  go.fisher <- runTest(go.topGO, algorithm="classic", statistic="fisher")
  go.ks <- runTest(go.topGO, algorithm="classic", statistic="ks")
  go.elim.ks <- runTest(go.topGO, algorithm="elim", statistic="ks")
  
  result <- GenTable(go.topGO, Fisher=go.fisher, KS = go.ks, elim.KS=go.elim.ks, topNodes=topNodes, orderBy=orderMethod)
  attr(result, "fisher") <- go.fisher
  attr(result, "ks") <- go.ks
  attr(result, "elim") <- go.elim.ks
  attr(result, "topGO") <- go.topGO
  
  #return result table
  return(result)
}