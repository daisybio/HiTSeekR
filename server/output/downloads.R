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
# output$downloadConsensusHits <- downloadHandler(
#   filename = function() { paste('consensus', 'hits', input$margin, input$method, '.csv', sep='_') },
#   content = function(file) {
#     data <- consensusHitList()
#     data$miRBase.url <- gsub("'", "", str_extract(data$Accession, "'http://.*?'"))
#     data$Accession <- gsub("<|>", "", str_extract(data$Accession, ">.*?<"))
#     data$signal <- NULL
#     data$signal.sem <- NULL
#     write.table(data, file, row.names=F, sep=",", quote=F)
#   }
# )

# mRNA targets #
output$downloadTargets <- downloadHandler(
  filename = function() { paste('mirna', 'target','genes', paste(input$selectedTargetDBs, collapse="_"), input$margin, input$method, input$normalization, '.csv', sep='_') },
  content = function(file) {
    data <- mirna.targets()
    write.table(data, file, row.names=F, sep=",", quote=F)
  }  
)

# drug targets #
output$downloadDrugTargets <- downloadHandler(
  filename = function() { paste('drug', 'target','genes', paste(input$selectedTargetDBs, collapse="_"), input$margin, input$method, input$normalization, '.csv', sep='_') },
  content = function(file) {
    data <- drug.targets()
    write.table(data, file, row.names=F, sep=",", quote=F)
  }  
)

# miRNA target gene permutation test result download #
output$downloadTargetPermutationTestResult <- downloadHandler(
  filename = function() { paste('mirna', 'target','genes', 'permutation', 'test', 'result', 
                                paste(input$selectedTargetDBs, collapse="_"), 
                                input$margin, 
                                input$method, 
                                input$normalization, '.csv', sep='_') 
                        },
  content = function(file) {
    data <- filtered.mirna.target.permutation()
    data <- data %>% dplyr::select(-Samples)
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
    if(input$screenType == "miRNA"){
      data <- targets.indicator.matrix()
    } 
    else if(input$screenType == "compound")
    {
      data <- drug.indicator.matrix()
    }
    else {
      data <- genes.indicator.matrix()
    } 
    write.table(data, file, sep="\t", quote=F, col.names=F)
  }
)

# Gene set analysis with htsanalyzer #
output$htsanalyzer.results.download.GO_CC <- downloadHandler(
  filename = function() { return(paste(input$htsanalyzer.resultType, "_GO_CC.txt", sep=""))},
  content = function(file) {
    data <- htsanalyzer.results()[["GO_CC"]]
    write.table(data, file, sep="\t", quote=F)
  }
) 
  
output$htsanalyzer.results.download.GO_MF <- downloadHandler(
  filename = function() { return(paste(input$htsanalyzer.resultType, "_GO_MF.txt", sep=""))},
  content = function(file) {
    data <- htsanalyzer.results()[["GO_MF"]]
    write.table(data, file, sep="\t", quote=F)
  }
  ) 

output$htsanalyzer.results.download.GO_BP <- downloadHandler(
  filename = function() { return(paste(input$htsanalyzer.resultType, "_GO_BP.txt", sep=""))},
  content = function(file) {
    data <- htsanalyzer.results()[["GO_BP"]]
    write.table(data, file, sep="\t", quote=F)
  }
  ) 

output$htsanalyzer.results.download.PW_KEGG <- downloadHandler(
  filename = function() { return(paste(input$htsanalyzer.resultType, "_KEGG.txt", sep=""))},
  content = function(file) {
    data <- htsanalyzer.results()[["PW_KEGG"]]
    write.table(data, file, sep="\t", quote=F)
  }
  ) 

output$htsanalyzer.results.download.REACTOME <- downloadHandler(
  filename = function() { return(paste(input$htsanalyzer.resultType, "_REACTOME.txt", sep=""))},
  content = function(file) {
    data <- htsanalyzer.results()[["REACTOME"]]
    write.table(data, file, sep="\t", quote=F)
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
