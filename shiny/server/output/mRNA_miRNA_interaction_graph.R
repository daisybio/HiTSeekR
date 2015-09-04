# plot a miRNA -> mRNA(genes) interaction graph
output$interactionGraph <- renderPlot({
  
  targets <- targetsForInteractionGraph()
  
  # variable that sets the least amount of miRNA's that it should connect.
  nrConnects <- 2
  
  # Using hashes for the lists, only requires average O(1) time for read/write
  
  #List of genes
  genes <- new.env( hash = TRUE, parent = emptyenv())
  
  #List of interactions including the gene
  geneInteraction <- new.env( hash = TRUE, parent = emptyenv())
  
  #List of category for miRNA
  miRNA_categories <- new.env( hash = TRUE, parent = emptyenv())
  
  print("-----------------")
  
  for(i in 1:nrow(targets)){
    gene_symbol <- as.character(targets[i, "gene_symbol"])
    miRNA <- as.character(targets[i, "miRNA"]) 
    
    # building hash for the number of observed connections, 
    # and the list of genes that it interacted with 
    if(is.null(genes[[gene_symbol]])){
      genes[[gene_symbol]] <- 1
      geneInteraction[[gene_symbol]] <- c(miRNA)
    }else{
      genes[[gene_symbol]] <- genes[[gene_symbol]] + 1     
      geneInteraction[[gene_symbol]] <- append(geneInteraction[[gene_symbol]] , miRNA)
    }
    
    # add the category to list of miRNA categories
    if(is.null(miRNA_categories[[miRNA]])){
      miRNA_categories[[miRNA]] <- targets[i, "category"]
    }
    
    NULL
  }
  
  # The edge list
  mat <- matrix(ncol = 2)
  
  
  # Build the edge list (matrix)
  for(v in ls(genes)){
    # for each gene found, see if the amount of connects is at least the minimum amount
    if(genes[[v]] >= nrConnects){
      for(gI in geneInteraction[[v]]){
        mat <- rbind(mat, c(v, gI))
      }
      
    }
  }
  
  # Must have at least 3 rows (+ 1 for [NA, NA] in matrix), else it fails
  if(nrow(mat) >= 4){ 
    colnames(mat) <- c("from", "to")
    mat = mat[-1,]
    print(mat)
    g <- qgraph(mat,directed=FALSE, layout = "spring")
    plot(g)
  }else{
    plot(1, type="n", axes=F, xlab="", ylab="")
    title(unlist(c("Too few genes with at least ", nrConnects, " connections were found.", "No graph will be drawn.")))
  }
  
})
