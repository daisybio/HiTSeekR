##This function computes observed and permutation-based scores associated 
##with a gene set enrichment analysis for a collection of Gene Sets.

my.collectionGsea <- function(collectionOfGeneSets, geneList, exponent=1, 
                              nPermutations=1000, minGeneSetSize=15, verbose=TRUE) {
  ##check input arguments
  paraCheck("gsc", collectionOfGeneSets)
  paraCheck("genelist", geneList)
  paraCheck("exponent", exponent)
  paraCheck("minGeneSetSize", minGeneSetSize)
  geneList.names <- names(geneList)
  paraCheck("nPermutations", nPermutations)  
  ##tag the gene sets that can be used in the analysis, i.e. those 
  ##that are smaller than the size of the gene list and that have more 
  ##than 'minGeneSetSize' elements that can be found in the geneList  
  nGeneSets <- length(collectionOfGeneSets)
  tagGeneSets <- rep(FALSE, nGeneSets)
  tagGeneSets[which(unlist(lapply(collectionOfGeneSets, length)) < 
                      length(geneList))] <- TRUE
  tagGeneSets[which(unlist(lapply(lapply(collectionOfGeneSets, 
                                         intersect, y=geneList.names), length)) < minGeneSetSize)] <- FALSE
  ##check that there are actually some gene sets that pass the max 
  ##and min cutoffs
  n.tagGeneSets <- sum(tagGeneSets)
  if(n.tagGeneSets == 0) 
    warning(paste("There are no gene sets in your collection",
                  " that pass the cutoffs on size", sep=""))
  if(n.tagGeneSets > 0) {
    ##Generate a matrix to store the permutation-based scores, with 
    ##one row for each gene set (that has been tagged) and one column 
    ##for each permutation	
    scoresperm <- matrix(rep(0, (nPermutations * n.tagGeneSets)), 
                         nrow=n.tagGeneSets)
    rownames(scoresperm) <- names(collectionOfGeneSets)[which(tagGeneSets)]
    ##Generate a vector to store the experimental scores
    ##one entry for each gene set (that has been tagged)
    scoresObserved <- rep(0, n.tagGeneSets)
    names(scoresObserved) <- names(collectionOfGeneSets)[which(tagGeneSets)]
    ##Compute the scores	
    ##create permutation gene list
    perm.gL <- sapply(1:nPermutations, function(n) names(geneList)[
      sample(1:length(geneList), length(geneList),replace=FALSE)])
    perm.gL<-cbind(names(geneList),perm.gL)
    ##check if package snow has been loaded and a cluster object 
    ##has been created for HTSanalyzeR	
    if(is(getOption("cluster"), "cluster") && 
         "package:snow" %in% search()) {
      scores <- gseaScoresBatchParallel(geneList, geneNames.perm = perm.gL,
                                        collectionOfGeneSets=collectionOfGeneSets[which(tagGeneSets)],
                                        exponent=exponent,nPermutations=nPermutations)
      sapply(1:n.tagGeneSets, function(i) {
        scoresperm[i,]<<-unlist(scores["scoresperm",i])
        scoresObserved[i]<<-unlist(scores["scoresObserved",i])
      }			
      )
    } else if(c("package:foreach") %in% search()){
      scores <- gseaScoresBatchForeach(geneList, geneNames.perm = perm.gL,
                                       collectionOfGeneSets=collectionOfGeneSets[which(tagGeneSets)],
                                       exponent=exponent,nPermutations=nPermutations)
      sapply(1:n.tagGeneSets, function(i) {
        scoresperm[i,]<<-unlist(scores[[i]][["scoresperm"]])
        scoresObserved[i]<<-unlist(scores[[i]][["scoresObserved"]])
      }  		
      )
    } else {
      if(verbose) 
        pb <- txtProgressBar(style=3)
      for(i in 1:n.tagGeneSets) {
        scores <- gseaScoresBatch(geneList, geneNames.perm=perm.gL, 
                                  geneSet=collectionOfGeneSets[[which(tagGeneSets)[i]]],
                                  exponent=exponent, nPermutations=nPermutations)
        scoresObserved[i] <- scores$scoresObserved
        scoresperm[i,] <- scores$scoresperm
        if(verbose) 
          setTxtProgressBar(pb, i/n.tagGeneSets)
      }	
      if(verbose) 
        close(pb)
    }
  } else {
    scoresObserved <- NULL
    scoresperm <- NULL
  }
  return(list("Observed.scores" = scoresObserved , "Permutation.scores" = scoresperm))	
}

R.utils::reassignInPackage("collectionGsea", "HTSanalyzeR", my.collectionGsea)
