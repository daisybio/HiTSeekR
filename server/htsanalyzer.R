htsanalyzer <- reactive({
  if(is.null(input$startHTSanalyzer)) return(NULL)
  if(input$startHTSanalyzer == 0) return(NULL)  
  
  #prepare progress bar
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = "Preparing input data")    
  
  isolate({      
      if(input$screenType %in% c("miRNA", "compound")) doGSEA <- FALSE
      else doGSEA <- input$htsanalyzer.doGSEA
      
      if(doGSEA && input$htsanalyzer.permutations > gsea.max.permutations)
      {
        showshinyalert(session, "htsanalyzer_status", paste("The current configuration allows a maximum of", gsea.max.permutations, "GSEA permutations. Please select a smaller value."), "danger")
        return(NULL)
      }
      if(input$screenType == "miRNA")              
      {        
          #get potential target genes from RNAhybrid. remove first entry (NA)
          if(input$selectedTargetDBs == "RNAhybrid_hsa"){
            pot.targets <- getRNAhybridTargetCounts(pval.threshold = input$rnah.p.value.threshold)
            gene.ids <- collect(pot.targets)$gene[-1]
          }
          else if(grepl("DIANA", input$selectedTargetDBs)){
            showshinyalert(session, "htsanalyzer_status", paste("Gene set analysis is not supported for DIANA miRNA targets since we do not know the universe size (number of potential targets)."), "danger")
            return(NULL)
          }
          else{
            pot.targets <- rmir.counts[[input$selectedTargetDBs]]
            gene.ids <- pot.targets$gene_id
          }
          all.samples.vector <- rep(1, length(gene.ids))
          names(all.samples.vector) <- gene.ids          
        if(input$htsanalyzer.miRNA.list == "miRNA_targets"){
          hit.list <- mirna.targets()
          hit.list <- as.character(hit.list$gene_id)
        }else if(input$htsanalyzer.miRNA.list == "miRNA_permutation")
        {
          if(input$mirna.target.permutation.button == 0){
            showshinyalert(session, "htsanalyzer_status", "You need to compute miRNA target gene significance first.", "danger")
            return(NULL)
          } 
          hit.list <- filtered.mirna.target.permutation()
          hit.list <- as.character(hit.list$gene_id)
        } else if(input$htsanalyzer.miRNA.list == "miRNA_KPM")
        {
          if(is.null(KPM.result())) {
            showshinyalert(session, "htsanalyzer_status", "You need to perform a KeyPathwayMiner enrichment analysis first.", "danger")
            return(NULL)
          }
          hit.list <- KPM.result()
          hit.list <- as.character(hit.list$gene.id) 
        }        
      }
      else if(input$screenType == "compound")
      {
         hit.list <- drug.targets()
         hit.list <- as.character(hit.list$gene_id)
         drug.targets <- drug.target.genes()
         all.samples.vector <- rep(1, length(drug.targets))
         names(all.samples.vector) <- drug.targets
      }
      else{        
        all.samples <- data()
        all.samples <- all.samples %>% dplyr::filter(Experiment %in% !!input$experimentSelected, Readout %in% !!input$readoutSelected)
                
        hit.list <- outliers()
        hit.list <- na.omit(hit.list$gene_id)
        
        all.samples.vector <- all.samples[,input$normalization]
        names(all.samples.vector) <- all.samples$gene_id
        #all.samples.vector <- all.samples.vector[-which(is.na(all.samples$gene_id))]
      }    
      
      ListGSC <- geneSets()              
      hideshinyalert(session, "htsanalyzer_status")  
      progress$set(message = "Preparing HTSanalyzeR for input data...")         
      gsca <- new("GSCA", listOfGeneSetCollections=ListGSC, geneList=all.samples.vector, hits=hit.list)
      gsca <- preprocess(gsca, species="Hs", initialIDs="Entrez.gene", keepMultipleMappings=TRUE, duplicateRemoverMethod="average", orderAbsValue=FALSE, progress=progress)
      gsca <- analyze(gsca, para=list(pValueCutoff=1, 
                                      pAdjustMethod = "BH",#input$htsanalyzer.adjust.method, 
                                      nPermutations=input$htsanalyzer.permutations, 
                                      minGeneSetSize=input$htsanalyzer.minimum.gene.set.size, 
                                      exponent=1), doGSEA=doGSEA, verbose=TRUE, progress=progress)
      return(gsca)
    })
})

htsanalyzer.annotated <- reactive({
  if(is.null(htsanalyzer())) return(NULL)
  geneset.types <- input$htsanalyzer.geneset.types
  ontologies <- setdiff(geneset.types, c("PW_KEGG", "REACTOME"))
  if(length(ontologies) == 0) ontologies <- NULL
  if("PW_KEGG" %in% geneset.types) kegg <- c("PW_KEGG")
  else kegg <- NULL
  
  result <- appendGSTerms(htsanalyzer(), goGSCs=ontologies, keggGSCs=kegg)  
  if("REACTOME" %in% geneset.types)
  {
    #extract from geneset collection where is was saved as attribute. remove species separated by colon
    path_names <- attr(result@listOfGeneSetCollections$REACTOME, "pathway_names")
    pathway_names <- str_split_fixed(path_names$path_name, ": ", 2)[,2]
    names(pathway_names) <- path_names$DB_ID
    if(!is.null(result@result$HyperGeo.results$REACTOME))
    {
      result@result$HyperGeo.results$REACTOME$Gene.Set.Term <- pathway_names[row.names(result@result$HyperGeo.results$REACTOME)]
    }
    if(!is.null(result@result$GSEA.results$REACTOME))
    {
      result@result$GSEA.results$REACTOME$Gene.Set.Term <- pathway_names[row.names(result@result$GSEA.results$REACTOME)]
    }
    if(!is.null(result@result$Sig.adj.pvals.in.both$REACTOME))
    {
      result@result$Sig.adj.pvals.in.both$REACTOME$Gene.Set.Term <- pathway_names[row.names(result@result$Sig.adj.pvals.in.both$REACTOME)]
    }
  }
  
  return(result)
})

htsanalyzer.results <- reactive({
  if(is.null(htsanalyzer.annotated())) return(NULL)
  gsca.annotated <- htsanalyzer.annotated()
  result <- gsca.annotated@result
  result <- result[[input$htsanalyzer.resultType]]
  lapply(result, function(x){
    if(input$htsanalyzer.resultType == "Sig.adj.pvals.in.both")
    {      
      dplyr::filter(x, HyperGeo.Adj.Pvalue <= input$htsanalyzer.pval.cutoff, GSEA.Adj.Pvalue <= input$htsanalyzer.pval.cutoff)     
    } else{
      dplyr::filter(x, Adjusted.Pvalue <= input$htsanalyzer.pval.cutoff)      
    }    
  })  
})

geneSets <- reactive({
  #species <- input$htsanalyzer.species #hardcoded for now
  species <- "Hs"
  geneset.types <- input$htsanalyzer.geneset.types
  ontologies <- str_replace_all(string = setdiff(geneset.types, c("PW_KEGG", "REACTOME")), pattern = "GO_", replacement = "")
  ontologies <- foreach(ontology = ontologies) %do% GOGeneSets(species = species, ontologies=ontology)
  if("PW_KEGG" %in% geneset.types) ontologies[[length(ontologies) + 1]] <- KeggGeneSets(species=species)
  if("REACTOME" %in% geneset.types) ontologies[[length(ontologies) + 1]] <- getReactomeGeneSets()
  names(ontologies) <- geneset.types
  #list(GO_MF=GO_MF, GO_BP=GO_BP, GO_CC=GO_CC, PW_KEGG=PW_KEGG)
  return(ontologies)
})