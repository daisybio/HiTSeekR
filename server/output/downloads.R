# Hit list #
output$downloadHits <- downloadHandler(
  filename = function() { paste('hits', input$normalization, input$margin, input$method, '.csv', sep='_') },
  content = function(file) {
    data <- outliers()
    data$miRBase.url <- gsub("'", "", str_extract(data$Accession, "'http://.*?'"))
    data$Accession <- gsub("<|>", "", str_extract(data$Accession, ">.*?<"))
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
    data$miRBase.url <- gsub("'", "", str_extract(data$Accession, "'http://.*?'"))
    data$Accession <- gsub("<|>", "", str_extract(data$Accession, ">.*?<"))
    data$signal <- NULL
    data$signal.sem <- NULL
    write.table(data, file, row.names=F, sep=",", quote=F)
  }
)

# mRNA targets #
output$downloadTargets <- downloadHandler(
  filename = function() { paste('mirna', 'target','genes', paste(input$selectedTargetDBs, collapse="_"), input$margin, input$method, input$normalization, '.csv', sep='_') },
  content = function(file) {
    data <- mirna.targets()
    write.table(data, file, row.names=F, sep=",", quote=F)
  }  
)

#hotnet2 heat file
output$downloadHotnetGeneList <- downloadHandler(
  filename = function() { paste('target', 'genes', 'heatscores', paste(input$selectedTargetDBs, collapse="_"), input$margin, input$method, input$normalization, '.txt', sep='_') },
  content = function(file) {
    data <- mirna.targets()
    data <- data %>% group_by(gene_symbol) %>% summarize(pvalue=sum(pvalue)) %>% mutate(pvalue=-log10(pvalue))
    write.table(data, file, col.names = F, row.names = F, sep="\t", quote = F)
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
