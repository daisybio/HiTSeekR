## helper functions ##

generateRandomHitList <- function(exp.data, numberOfHits){  
  return(exp.data %>% sample_n(numberOfHits))  
}

generateIndicatorMatrix <- function(mirna.targets){
  indicator.matrix <- as.data.frame.matrix(table(mirna.targets[,c("gene_id", "mature_miRNA")]))
  indicator.matrix[indicator.matrix > 1] <- 1
  return(indicator.matrix)  
}

rbind.with.progress <- function(progress, steps){  
    count <- 0
    function(...) {
      new.entries <- ((length(list(...)) - 1) / steps) 
      count <<- count + new.entries
      progress$set(count, detail = paste(round(count * 100,2), "% done"))            
      rbind(...)
    }
}

# selectedHitList <- reactive({
#   if(input$useConsensus == "hit list") data <- outliers()
#   else data <- consensusHitList()
#   return(data)
# })

rnah.p.value.threshold <- reactive({
  if(is.null(input$rnah.p.value.threshold)) rnah.p.value.threshold <- NULL
  else rnah.p.value.threshold <- input$rnah.p.value.threshold
  return(rnah.p.value.threshold)
})

#compute % of how many miRNA members of any represented family are hit. 
family.hitrate <- reactive({ 
  all.data <- data()
  all.data <- all.data %>% dplyr::group_by(prefam_acc) %>% dplyr::summarise(library_count=dplyr::n_distinct(mirna_id)) 
  hits <- outliers()
  
  result <- formattedTable(hits, FALSE)
  result <- within(result, hits <- paste(category, Sample))
  result <- result %>% dplyr::group_by(prefam_acc, id) %>% dplyr::summarise(hits_count=n(),
                                                                            distinct_hits_count=n_distinct(mirna_id), 
                                                                            samples=paste(hits, collapse="")
  )
  
  final_result <- dplyr::left_join(result, all.data, by="prefam_acc")
  final_result <- final_result %>% mutate(family_coverage=distinct_hits_count/library_count)
  final_result <- final_result %>% filter(library_count > input$family_size_cutoff, family_coverage > (input$family_coverage_cutoff/100))
  if(nrow(final_result) == 0) stop("no families found")
  final_result <- as.data.frame(final_result)
  family_coverage <- as.numeric(final_result$family_coverage)
  
  #add color code for <0.33, 0.33-0.66, >0.66.
  final_result$family_coverage <- paste("<div style='background:#FDB462; text-align:center; border-radius: 15px; width:50px; height:25px;'>", 
                                        round(100*family_coverage, 0), "%",
                                        "</div>", sep="")
  final_result[family_coverage < 0.33, "family_coverage"] <-  sub("#FDB462", "#FF0000", final_result[family_coverage < 0.33, "family_coverage"])
  final_result[family_coverage > 0.66, "family_coverage"] <-  sub("#FDB462", "#40FF00", final_result[family_coverage > 0.66, "family_coverage"])
  final_result <- na.omit(final_result)
  
  return(as.data.frame(final_result))
})


mirna.hits <- reactive({
  data <- outliers()  
  if(is.null(data)) return(NULL)
  #result <- isolate({
  
  #check if accession type is mirna_id we need to convert to mature mirna ids in that case
  if(input$accessionColType == "MI")
  {
    showshinyalert(session, "mirna_target_status", "For pri-miRNA (Accession MI) it is unclear which of the two strands (3p or 5p) is active. Therefore, targets for both strands are predicted.","warning")
    mimat <- as.data.frame(mirbaseMATURE)
    data <- left_join(data, mimat, by=c("mirna_id"))
  }
  return(data)
})
## find miRNA target genes ##
mirna.targets <- reactive({
  
 data <- mirna.hits()
 
 #prepare progress bar
 progress <- shiny::Progress$new()
 on.exit(progress$close())
 progress$set(message = "Querying for miRNA target genes")    
 
  gene.symbols <- FALSE
  
  if(input$selectedTargetDBs != "RNAhybrid_hsa"){
      gene.symbols <- TRUE    
  }
 
  tryCatch({ 
    result <- getTargets(data, rnah.pvalue.threshold=rnah.p.value.threshold(), get.gene.symbols = gene.symbols,databases=input$selectedTargetDBs, diana.threshold = input$diana.microT.min.score)
  }, error = function(e){ 
      showshinyalert(session, "mirna_target_status", paste("An error occured: ", e$message, sep="") ,"danger")
      return(NULL)
    })
  
  if(nrow(result) == 0) stop("No target has been found. reduce stringency or increase number of hits.")
  
  return(result)

})

rnah.mirna.count <- reactive({
  getRNAhybridNumOfmiRNAs()
})

observe({
  if(is.null(input$mirna.target.permutation.button)) return(NULL)
  if(input$selectedTargetDBs != "RNAhybrid_hsa" && input$highConfidenceTargetsMethod == "hypergeometric test")
  {
    showshinyalert(session, "mirna_conf_status", "Hypergeometric test is only supported for RNAhybrid.","warning")
  }
  else if(grepl("DIANA", input$selectedTargetDBs, ignore.case=TRUE)){
    showshinyalert(session, "mirna_conf_status", "High-confidence miRNA target genes can not be computed with webservices due to excessive querying during the permutation test.","warning")
  }
  else if(input$mirna.target.permutations > gsea.max.permutations){
    showshinyalert(session, "mirna_conf_status", paste("The current configuration of HiTSeekR only allows a maximum of", gsea.max.permutations, " permutations. The selected value will be adjusted accordingly.", sep = ""),"danger")
  }
  else{
    hideshinyalert(session, "mirna_conf_status")
  }
  
})

mirna.target.permutation <- eventReactive(input$mirna.target.permutation.button,{

    if(input$selectedTargetDBs != "RNAhybrid_hsa" && input$highConfidenceTargetsMethod == "hypergeometric test")
    {
      return(NULL)
    }
    else if(grepl("DIANA", input$selectedTargetDBs, ignore.case=TRUE)){
      return(NULL)
    }
    if(input$mirna.target.permutations > gsea.max.permutations){
      permutations <- gsea.max.permutations
    }
    else{
      permutations <- input$mirna.target.permutations
    }
    
    #get hit list and mirna targets
    hit.list <- outliers()
    hit.list <- formattedTable(hit.list, FALSE)
    hit.list <- within(hit.list, Sample <- paste(category, Sample))
    
    if(input$accessionColType == "MI")
    {
      mimat <- as.data.frame(mirbaseMATURE)
      hit.list <- left_join(hit.list, mimat, by=c("mirna_id"))
    }
    
    hit.list.targets <- as.data.frame(mirna.targets())
    
    #combine original hit list with original target list to get number of distinct mirna families and miRNA names    
    hit.list.targets.combined <- left_join(hit.list.targets, hit.list, by=c("mature_miRNA"= "mature_name"))  
    target.list.by.gene <- hit.list.targets.combined %>% group_by(gene_id, gene_symbol) %>% dplyr::summarise(number_of_miRNAs = n_distinct(mature_miRNA), mirna_families = n_distinct(prefam_acc), Samples=paste(Sample, collapse="\n"), mature_miRNA=paste(mature_miRNA, collapse = "/"))  
    
    #focus only on genes that are targets of the original hit list
    genes.of.interest <- unique(hit.list.targets$gene_id)
    
    #prepare progress bar
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    
    if(input$highConfidenceTargetsMethod == "permutation test")
    {
      exp.data <- data() %>% filter(!is.na(Accession))
      exp.data <- as.data.frame(exp.data) 
      rnah.p.value.threshold <- rnah.p.value.threshold()
      selectedTargetDBs <- input$selectedTargetDBs

      progress$set(message = "Generating random hit lists")    
      
      #randomization loop
      result <- foreach(i = 1:permutations, .combine=rbind.with.progress(progress, permutations),
                        .packages=c("dplyr", "RmiR", "foreach"), 
                        .export=c("generateRandomHitList", "getTargets", "data.folder", "getRNAhybridTargets")) %dopar%
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
      progress$set(message = "Summarizing results")    
      
      result <- result %>% dplyr::group_by(gene_id) %>% dplyr::summarise(n= sum(n)) %>% mutate(n = n / permutations) %>% dplyr::rename(expected = n)              
      
      #combine information from the processed hit list with the result
      mirna.hit.counts <- left_join(target.list.by.gene, result, by="gene_id") %>%
        dplyr::rename(count=expected)     
    }
    else{      
      progress$set(message = "Performing hypergeometric tests")    
      mirna.hit.counts <- getRNAhybridTargetCounts(genes.of.interest, input$rnah.p.value.threshold)      
      mirna.hit.counts <- left_join(target.list.by.gene, mirna.hit.counts, by=c("gene_id" = "gene"), copy = TRUE)
    }  
    #universe and number of hits
    universe.size <- 1242 #rnah.mirna.count()
    hit.list.size <- length(unique(hit.list.targets$mature_miRNA))
    
    #calculate cumulative hypergeometric test, e.g. P(X >= k). We use number_of_miRNAs - 1, because otherwise it would be P(X > k)
    final_result <- mirna.hit.counts %>% mutate(p.value = phyper(number_of_miRNAs-1, count, universe.size - count, hit.list.size, lower.tail=FALSE))          
    progress$set(message = "Ajusting p-values (BH)")
    final_result$p.adj <- p.adjust(final_result$p.value, method="BH")
                                            
    return(final_result)
  })

filtered.mirna.target.permutation <- reactive({
  result <- mirna.target.permutation()
  if(is.null(result)) return(NULL)
     
  final_result <- result %>% filter(number_of_miRNAs >= input$mirna.target.permutation.num.of.mirnas.cutoff) %>%
    filter(p.adj < input$mirna.target.permutation.padj.cutoff)
  return(final_result)
})

list.of.random.mirna.indicator.matrices <- reactive({
  if(input$accessionColType != "MIMAT") stop("This function is only supported for datasets providing MIMAT ids")

  rnah.p.value.threshold <- rnah.p.value.threshold()
  
  result <- foreach(i = 1:input$random.miRNA.iterations) %dopar% {    
    data <- generateRandomHitList(length(unique(outliers()$Accession)))
    mi.targets <- getTargets(data, rnah.pvalue.threshold=rnah.p.value.threshold, databases=input$selectedTargetDBs)  
    ind.matrix <- generateIndicatorMatrix(mi.targets)    
    return(ind.matrix)
  }
})

#indicator_matrix 
targets.indicator.matrix <- reactive({
  mi.targets <- mirna.targets() 
  if(input$KPM.miRNA.list != "miRNA_targets")
  { 
    if(input$mirna.target.permutation.button == 0){
      showshinyalert(session, "kpm_status", "You need to determine high confidence miRNA target genes first. Go to the microRNA tab.", "danger")   
      return(NULL)
    } 
    else if(nrow(filtered.mirna.target.permutation()) == 0){
      showshinyalert(session, "kpm_status", "No high confidence miRNA target genes were found with the selected settings. You need to change these settings before you can continue here.", "danger")   
      return(NULL)
    }
    confidence.genes <- filtered.mirna.target.permutation()
    mi.targets <- mi.targets %>% filter(mi.targets$gene_id %in% confidence.genes$gene_id)  
  }
  generateIndicatorMatrix(mi.targets)
})

ind.matrix.props <- renderText({
  paste(nrow(targets.indicator.matrix()), "genes\n", ncol(targets.indicator.matrix()), "cases")  
})

# DIANA mirpath results
mirpath.results <- reactive({
  
  #prepare progress bar
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = "Querying DIANA mirPATH. Waiting for results...")   
  
  get_DIANA_mirPath(miRNAs = mirna.hits(), threshold = input$mirpath_threshold, geneIntersectionCutoff = input$mirpath_cutoff, selection = input$mirpath_selection, fdr = input$mirpath_fdr,conservative = input$mirpath_conservative)
})