gseaScoresBatchForeach <- function (geneList, geneNames.perm, collectionOfGeneSets, exponent = 1, 
                                    nPermutations = 1000) 
{
  paraCheck <- HTSanalyzeR:::paraCheck
  paraCheck("genelist", geneList)
  paraCheck("gsc", collectionOfGeneSets)
  paraCheck("exponent", exponent)
  paraCheck("nPermutations", nPermutations)
  if (!is.matrix(geneNames.perm)) 
    stop("'geneNames.perm' should be a matrix!\n")
  if (ncol(geneNames.perm) != (nPermutations + 1)) 
    stop("The No of columns of 'geneNames.perm' should be equal to 'nPermutations'!\n")
  gseaScoresBatchLocal <- function(geneList, geneNames.perm, 
                                   geneSet, exponent, nPermutations) {
    geneList.names <- names(geneList)
    geneSet <- intersect(geneList.names, geneSet)
    nh <- length(geneSet)
    N <- length(geneList)
    ES <- rep(0, nPermutations + 1)
    Phit <- matrix(0, nrow = N, ncol = nPermutations + 1)
    Pmiss <- Phit
    runningES <- NULL
    if (nh > N) 
      stop("Gene Set is larger than Gene List")
    hits <- matrix(FALSE, nrow = N, ncol = nPermutations + 
                     1)
    hits[which(!is.na(match(geneNames.perm, geneSet)))] <- TRUE
    hits <- matrix(hits, ncol = nPermutations + 1, byrow = FALSE)
    if (sum(hits[, 1]) > 0) {
      junk <- sapply(1:(nPermutations + 1), function(i) Phit[which(hits[, 
                                                                        i]), i] <<- abs(geneList[which(hits[, i])])^exponent)
      NR <- colSums(Phit)
      Pmiss[which(!hits)] <- 1/(N - nh)
      Pmiss <- sapply(1:(nPermutations + 1), function(i) cumsum(Pmiss[, 
                                                                      i]))
      Phit <- sapply(1:(nPermutations + 1), function(i) cumsum(Phit[, 
                                                                    i])/NR[i])
      runningES <- Phit - Pmiss
      ESrange <- sapply(1:(nPermutations + 1), function(i) range(runningES[, 
                                                                           i]))
      ES <- sapply(1:(nPermutations + 1), function(i) ESrange[which.max(abs(ESrange[, 
                                                                                    i])), i])
      if (is.list(ES)) 
        ES <- unlist(ES)
    }
    ES <- list(scoresObserved = ES[1], scoresperm = ES[2:(nPermutations + 
                                                            1)])
    return(ES)
  }
  scores <- foreach(i=1:length(collectionOfGeneSets), .export=c("paraCheck")) %dopar% 
{
  gseaScoresBatchLocal(geneList, geneNames.perm = geneNames.perm, 
                       geneSet = as.integer(collectionOfGeneSets[[i]]), 
                       exponent = exponent, nPermutations = nPermutations)
}
return(scores)
}