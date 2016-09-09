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

# miRNA targets #
output$downloadTargets <- downloadHandler(
  filename = function() { paste('mirna', 'target','genes', paste(input$selectedTargetDBs, collapse="_"), input$margin, input$method, input$normalization, '.txt', sep='_') },
  content = function(file) {
    data <- mirna.targets()
    write.table(data, file, row.names=F, col.names = T, sep="\t", quote=F)
  }  
)

#mirPath results #
output$downloadMirPathResults <- downloadHandler(
  filename = function() { "mirpath_results.txt" },
  content = function(file) {
    data <- mirpath.results()
    write.table(data, file, row.names=F, col.names = T, sep="\t", quote=F)
  }  
)

# drug targets #
output$downloadDrugTargets <- downloadHandler(
  filename = function() { paste('drug', 'target','genes', paste(input$selectedTargetDBs, collapse="_"), input$margin, input$method, input$normalization, '.txt', sep='_') },
  content = function(file) {
    data <- drug.targets()
    write.table(data, file, row.names=F, sep="\t", quote=F)
  }  
)

# miRNA target gene permutation test result download #
output$downloadTargetPermutationTestResult <- downloadHandler(
  filename = function() { paste('mirna', 'target','genes', 'permutation', 'test', 'result', 
                                paste(input$selectedTargetDBs, collapse="_"), 
                                input$margin, 
                                input$method, 
                                input$normalization, '.txt', sep='_') 
                        },
  content = function(file) {
    data <- filtered.mirna.target.permutation()
    data <- data %>% dplyr::select(-Samples)
    write.table(data, file, row.names=F, sep="\t", quote=F)
  }  
)

# mircancer results
output$downloadMirCancerDbResult <- downloadHandler(
  filename = function() { paste('mircancerdb', 'result', 
                                input$margin, 
                                input$method, 
                                input$normalization, '.txt', sep='_') 
  },
  content = function(file) {
    data <- hits.mircancer()
    write.table(data, file, row.names=F, sep="\t", quote=F)
  }  
)


# mir family results
output$downloadMirFamilyResult <- downloadHandler(
  filename = function() { paste('mirfamilies', 'result', 
                                input$margin, 
                                input$method, 
                                input$normalization, '.txt', sep='_') 
  },
  content = function(file) {
    data <- family.hitrate()
    
    data$samples <- str_replace_all(data$samples, "<div.+?>S</div> ", "(S)")
    data$samples <- str_replace_all(data$samples, "<div.+?>P</div> ", "(P)")
    data$samples <- str_replace_all(data$samples, "<div.+?>I</div> ", "(I)")
    
    write.table(data, file, row.names=F, sep="\t", quote=F)
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
    write.table(data, file, sep="\t", quote=F, col.names=T)
  }
)

# Download currently selected KPM result graph as SIF file
output$download_kpm_SIF <- downloadHandler(
  filename = function() { 
    if(input$kpm_union_graph) return(paste("kpm_result_union_graph.sif", sep=""))
    else return(paste("kpm_result_graph_", input$kpm_selected_solution, ".sif", sep=""))},
  content = function(file) {
    edges <- kpm.graph.data()[[1]]
    edges <- cbind(edges[,1], 'pp', edges[,2])
    if(input$screenType == "miRNA") {
      edges[-grep("^[0-9]+$", edges[,1]), 2] <- "mi"
    }
    if(input$screenType == "compound"){
      edges[-grep("^[0-9]+$", edges[,1]), 2] <- "di"
    }
    write.table(edges, file, sep="\t", row.names = FALSE, col.names = FALSE, quote=F)
  }
) 

# Download currently selected KPM result graph as table
output$download_kpm_node_table <- downloadHandler(
    filename = function(){
  if(input$kpm_union_graph) return(paste("kpm_result_union_graph.txt", sep=""))
  else return(paste("kpm_result_graph_", input$kpm_selected_solution, ".txt", sep=""))},
  content = function(file) {
    data <- kpm.node.table()
    write.table(data, file, sep="\t", row.names = FALSE, quote=F)
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
  filename = function() { return("GO.txt")},
  content = function(file) {
    data <- goEnrichment()
    write.table(data, file, row.names=F, sep="\t", quote=F)
  }
)

output$datasetName <- renderPrint({
  input$updateNormalization
  isolate(input$dataset)
})
