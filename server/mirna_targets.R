## helper functions ##

generateRandomHitList <- function(exp.data, numberOfHits){  
  return(exp.data %>% sample_n(numberOfHits))  
}

generateIndicatorMatrix <- function(mirna.targets){
  indicator.matrix <- as.data.frame.matrix(table(mirna.targets[,c("gene_id", "mature_miRNA")]))
  indicator.matrix[indicator.matrix > 1] <- 1
  return(indicator.matrix)  
}

rbind.with.progress <- function(progress, permutations){  
    function(...) {
      new.entries <- (length(list(...)) - 1) / permutations
      progress$inc(new.entries, detail = "please wait...")            
      rbind(...)
    }
}

rnah.p.value.threshold <- reactive({
  if(is.null(input$rnah.p.value.threshold)) rnah.p.value.threshold <- NULL
  else rnah.p.value.threshold <- 10^input$rnah.p.value.threshold
  return(rnah.p.value.threshold)
})

mirna.target.permutation.padj.cutoff <- reactive({
  if(is.null(input$mirna.target.permutation.padj.cutoff)) mirna.target.permutation.padj.cutoff <- NULL
  else mirna.target.permutation.padj.cutoff <- 10^input$mirna.target.permutation.padj.cutoff
  return(mirna.target.permutation.padj.cutoff)
})

## find miRNA target genes ##
mirna.targets <- reactive({
  
  input$updateTargets
  data <- selectedHitList()  
  
  result <- isolate({
    
    #check if accession type is mirna_id we need to convert to mature mirna ids in that case
    if(input$accessionColType == "MI")
    {
      showshinyalert(session, "mirna_target_status", "For pri-miRNA (Accession MI) it is unclear which of the two strands (3p or 5p) is active. Therefore, targets for both strands are predicted.","warning")
      mimat <- as.data.frame(mirbaseMATURE)
      data <- left_join(data, mimat, by=c("mirna_id"))
    }
    
    result <- getTargets(data, rnah.pvalue.threshold=rnah.p.value.threshold(), databases=input$selectedTargetDBs)
    #if(input$excludeDBcol) result <- result[,setdiff(colnames(result), c("db_list"))]
    
    if(nrow(result) == 0) stop("No target has been found. reduce stringency or increase number of hits.")
    
    return(result)
  })                 
  return(result)
})

mirna.target.permutation <- reactive({
  
  if(input$mirna.target.permutation.button == 0) return(NULL)
  
  isolate({
    #get hit list and mirna targets
    hit.list <- selectedHitList()
    hit.list.targets <- as.data.frame(mirna.targets())
    exp.data <- data() %>% filter(!is.na(Accession))
    exp.data <- as.data.frame(exp.data) 
    rnah.p.value.threshold <- rnah.p.value.threshold()
    selectedTargetDBs <- input$selectedTargetDBs
    permutations <- input$mirna.target.permutations
    
    #focus only on genes that are targets of the original hit list
    genes.of.interest <- unique(hit.list.targets$gene_id)
    
    #prepare progress bar
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "miRNA target permutation test")    
    
    #permutation loop
    result <- foreach(i = 1:permutations, .combine=rbind.with.progress(progress, permutations),
                      .packages=c("dplyr", "RmiR", "foreach"), 
                      .export=c("generateRandomHitList", "getTargets")) %dopar%
    {
      #get a random hit list of equal length
      test.data <- generateRandomHitList(exp.data, length(unique(hit.list$Accession)))    
      
      #get targets of the random hit list
      mi.targets <- getTargets(test.data, rnah.pvalue.threshold=rnah.p.value.threshold, databases=selectedTargetDBs)      
      mi.targets <- mi.targets %>% filter(gene_id %in% genes.of.interest)
      
      #count for each gene how often it was included
      result <- mi.targets %>% count(gene_id)    
      return(result)
    }
    #count for all permutations how often a gene was a target,
    #divide by the number of permutations to get an E-value
    
    result <- result %>% group_by(gene_id) %>% summarise(n= sum(n)) %>% mutate(n = n / permutations) %>% dplyr::rename(expected = n)        
    
    #combine original hit list with original target list to get number of distinct mirna families and miRNA names
    hit.list.targets.combined <- left_join(hit.list.targets, hit.list, by=c("mature_miRNA"= "mature_name"))  
    target.list.by.gene <- hit.list.targets.combined %>% group_by(gene_id, gene_symbol) %>% summarise(number_of_miRNAs = n_distinct(mature_miRNA), mirna_families = n_distinct(prefam_acc), mature_miRNA=paste(mature_miRNA, collapse = "/"))  
    
    #combine information from the step above with the result
    final_result <- left_join(target.list.by.gene, result, by="gene_id")
    final_result <- dplyr::select(final_result, gene_id, gene_symbol, number_of_miRNAs, expected_number_of_miRNAs=expected, mirna_families)        
    
    #calculate p-values using the Poisson distribution (assuming independence)
    final_result <- final_result %>% mutate(p.value = ppois(number_of_miRNAs, lambda=expected_number_of_miRNAs, lower=FALSE))
    final_result$p.adj <- p.adjust(final_result$p.value)
                                            
    return(final_result)
  })
})

filtered.mirna.target.permutation <- reactive({
  result <- mirna.target.permutation()
  if(is.null(result)) return(NULL)
     
  final_result <- result %>% filter(number_of_miRNAs > input$mirna.target.permutation.num.of.mirnas.cutoff) %>%
    filter(p.adj < mirna.target.permutation.padj.cutoff())
  return(final_result)
})

list.of.random.mirna.indicator.matrices <- reactive({
  if(input$accessionColType != "MIMAT") stop("This function is only supported for datasets providing MIMAT ids")

  if(is.null(input$rnah.p.value.threshold)) rnah.p.value.threshold <- NULL
  else rnah.p.value.threshold <- 10^input$rnah.p.value.threshold
  
  result <- foreach(i = 1:input$random.miRNA.iterations) %dopar% {    
    data <- generateRandomHitList(length(unique(selectedHitList()$Accession)))
    mi.targets <- getTargets(data, rnah.pvalue.threshold=rnah.p.value.threshold, databases=input$selectedTargetDBs)  
    ind.matrix <- generateIndicatorMatrix(mi.targets)    
    return(ind.matrix)
  }
})

#indicator_matrix 
targets.indicator.matrix <- reactive({
  mi.targets <- mirna.targets()
  generateIndicatorMatrix(mi.targets)
})