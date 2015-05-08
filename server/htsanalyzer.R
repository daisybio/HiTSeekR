htsanalyzer <- reactive({
  if(is.null(input$startHTSanalyzer)) return(NULL)
  if(input$startHTSanalyzer == 0) return(NULL)
  
  #prepare progress bar
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = "Preparing input data")    
  
  isolate({      
      if(input$screenType == "miRNA")              
      {                     
          #get potential target genes from RNAhybrid. remove first entry (NA)
          pot.targets <- getRNAhybridTargetCounts()
          
          gene.ids <- collect(pot.targets)$gene[-1]
          all.samples.vector <- rep(1, length(gene.ids))
          names(all.samples.vector) <- gene.ids          
        if(input$htsanalyzer.miRNA.list == "miRNA_targets"){
          hit.list <- mirna.targets()
          hit.list <- as.character(hit.list$gene_id)
        }else if(input$htsanalyzer.miRNA.list == "miRNA_permutation")
        {
          if(is.null(filtered.mirna.target.permutation())) stop("You need to perform a miRNA target permutation test first.")
          hit.list <- filtered.mirna.target.permutation()
          hit.list <- as.character(hit.list$gene_id)
        } else if(input$htsanalyzer.miRNA.list == "miRNA_KPM")
        {
          if(is.null(KPM.result())) stop("You need to perform a KeyPathwayMiner enrichment analysis first.")
          hit.list <- KPM.result()
          hit.list <- as.character(hit.list$gene.id) 
        }        
      } else{
        all.samples <- data()
        hit.list <- htsanalyzerHitList()
        hit.list <- na.omit(hit.list$gene_id)
        
        all.samples.vector <- all.samples[,input$normalization]
        names(all.samples.vector) <- all.samples$gene_id
        all.samples.vector <- all.samples.vector[-which(is.na(all.samples$gene_id))]
      }    
      
      ListGSC <- geneSets()
      
      if(input$screenType == "miRNA") doGSEA <- FALSE
      else doGSEA <- input$htsanalyzer.doGSEA
      progress$set(message = "Performing tests...")    
      gsca <- new("GSCA", listOfGeneSetCollections=ListGSC, geneList=all.samples.vector, hits=hit.list)
      gsca <- preprocess(gsca, species=input$htsanalyzer.species, initialIDs="Entrez.gene", keepMultipleMappings=TRUE, duplicateRemoverMethod="average", orderAbsValue=FALSE)
      gsca <- analyze(gsca, para=list(pValueCutoff=input$htsanalyzer.pval.cutoff, 
                                      pAdjustMethod =input$htsanalyzer.adjust.method, 
                                      nPermutations=input$htsanalyzer.permutations, 
                                      minGeneSetSize=input$htsanalyzer.minimum.gene.set.size, 
                                      exponent=1), doGSEA=doGSEA)
      return(gsca)
    })
})

htsanalyzer.annotated <- reactive({
  if(is.null(htsanalyzer())) return(NULL)
  geneset.types <- input$htsanalyzer.geneset.types
  ontologies <- setdiff(geneset.types, "PW_KEGG")
  if(length(ontologies) == 0) ontologies <- NULL
  if("PW_KEGG" %in% geneset.types) kegg <- c("PW_KEGG")
  else kegg <- NULL
  appendGSTerms(htsanalyzer(), goGSCs=ontologies, keggGSCs=kegg)  
})

htsanalyzer.results <- reactive({
  if(is.null(htsanalyzer.annotated())) return(NULL)
  gsca.annotated <- htsanalyzer.annotated()
  result <- gsca.annotated@result
  result <- result[[input$htsanalyzer.resultType]]
  return(result)
})

geneSets <- reactive({
  species <- input$htsanalyzer.species 
  geneset.types <- input$htsanalyzer.geneset.types
  ontologies <- str_replace_all(string = setdiff(geneset.types, "PW_KEGG"), pattern = "GO_", replacement = "")
  ontologies <- foreach(ontology = ontologies) %do% GOGeneSets(species = species, ontologies=ontology)
  if("PW_KEGG" %in% geneset.types) ontologies[[length(geneset.types)]] <- KeggGeneSets(species=species)
  names(ontologies) <- geneset.types
  #list(GO_MF=GO_MF, GO_BP=GO_BP, GO_CC=GO_CC, PW_KEGG=PW_KEGG)
  return(ontologies)
})