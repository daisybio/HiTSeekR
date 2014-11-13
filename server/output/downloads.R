# Hit list #
output$downloadHits <- downloadHandler(
  filename = function() { paste('hits', input$normalization, input$margin, input$method, '.csv', sep='_') },
  content = function(file) {
    data <- outliers()
    data$miRBase.url <- gsub("'", "", str_extract(data$miRBase.accession, "'http://.*?'"))
    data$miRBase.accession <- gsub("<|>", "", str_extract(data$miRBase.accession, ">.*?<"))
    data$signal <- NULL
    data$signal.sem <- NULL
    write.table(data, file, row.names=F, sep=",", quote=F)
  }
)

# Consensus hit list #
output$downloadConsensusHits <- downloadHandler(
  filename = function() { paste('consensus', 'hits', input$margin, input$method, '.csv', sep='_') },
  content = function(file) {
    data <- consensusHitList()
    data$miRBase.url <- gsub("'", "", str_extract(data$miRBase.accession, "'http://.*?'"))
    data$miRBase.accession <- gsub("<|>", "", str_extract(data$miRBase.accession, ">.*?<"))
    data$signal <- NULL
    data$signal.sem <- NULL
    write.table(data, file, row.names=F, sep=",", quote=F)
  }
)

# mRNA targets #
output$downloadTargets <- downloadHandler(
  filename = function() { paste('targets', paste(input$selectedTargetDBs, collapse="_"), input$margin, input$method, '.csv', sep='_') },
  content = function(file) {
    data <- targets()
    if(input$colorizeInTargetList){
      data$categories <- gsub("blue", "P", gsub("red", "S", sapply(str_extract_all(data$miRNA_list, "red|blue"), paste, collapse="/")))
      data$miRNA_list <- gsub("<|>", "", sapply(str_extract_all(data$miRNA_list, ">.*?<"), paste, collapse="/"))
    } 
    data$url <- str_extract(data$gene_symbol, "http://.*?%5D")
    data$gene_symbol <- gsub("<|>", "", str_extract(data$gene_symbol, ">.*?<"))
    
    write.table(data, file, row.names=F, sep=",", quote=F)
  }  
)

#indicator matrix
output$downloadIndicatorMatrix <- downloadHandler(
  filename = function() { return("indicator_matrix.txt")},
  content = function(file) {
    data <- targets.indicator.matrix()
    write.table(data, file, sep="\t", quote=F)
  }
)

# GO enrichment graph #
output$dlGnOntGraph <- downloadHandler(
  filename = function() { return("GO.pdf")},
  content = function(file) {
    data <- goEnrichment()
    pdf(file)
    showSigOfNodes(attr(data, "topGO"), score(attr(data, input$goSelectedMethod)), firstSigNodes = input$goSelectedNodes, useInfo=input$goUseInfo)
    dev.off()
  }
)

# GO enrichment table #
output$dlGnOntTbl <- downloadHandler(
  filename = function() { return("GO.csv")},
  content = function(file) {
    data <- goEnrichment()
    write.table(data, file, row.names=F, sep=",", quote=F)
  }
)

output$datasetName <- renderPrint({
  input$updateNormalization
  isolate(input$dataset)
})
