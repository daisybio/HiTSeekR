bayesianHitSelection <-function(dataTable, neg.ctrl="NEG", signalColumn="Raw", alpha=0.05, updateProgress=NULL){
  
  library(foreach)
  library(iterators)
  
  if(is.function(updateProgress))
    updateProgress(detail="Bayesian statistics", value=0)
  
  #Log transform values if dealing with raw values
  if(signalColumn=="Raw") values <- log2(dataTable$Raw)
  else values <- dataTable[[signalColumn]]
  X = data.frame(plate=as.integer(dataTable$Plate),row=dataTable$Row, column=dataTable$Column, value=values, control=dataTable$Control)
  
  #nr of plates
  A = max(X$plate)
  
  #transformed values of negative control plates
  Yij <- subset(X, control %in% neg.ctrl)
  
  #platewise mean negative control values
  Ydotj <- c()
  for(i in 1:A){
    ithSubset <- subset(X, plate==i & control %in% neg.ctrl)
    
    if(nrow(ithSubset) > 0){
      Ydotj[i] <- mean(ithSubset$value)
    }else{
      if(i > 1) Ydotj[i] <- Ydotj[i-1]
      else stop("Cannot set plate negative control mean")
    }
  }
  
  #difference of value from mean
  dij <- c()
  for(rowCount in 1:nrow(Yij)){
    currentY <- Yij[rowCount,]
    j <- currentY$plate
    dij[rowCount] <- currentY$value - Ydotj[j]
  }
  
  #estimate of variance of negative controls (nrow(Yij) == sum of control values)
  sigmaSqrd  <- sum(dij^2)/(nrow(Yij)-A)
  
  tauSqrd <- max(var(X$value, na.rm=T)-sigmaSqrd,0)
  
  #estimated variance of all Xi
  varX <- sigmaSqrd + tauSqrd
  rowCount <- 1
  rowTotal <- nrow(X)
  rowStep <- floor(rowTotal / 10)
  rowNextStep <- rowStep
  percentage <- 0
  
  pHs <- foreach(row=iter(X, by='row'), .combine=rbind.with.progress(updateProgress, rowTotal), .verbose=FALSE, .export = c("bayesianHypothesisTesting")) %dopar%
  {     
    #for each sample we need a prior distribution
    prior_mean <- Ydotj[row$plate]
    prior_var <- tauSqrd    
    
    #In the posterior distribution we can test our three hypothesis
    posterior_mean <- (Ydotj[row$plate]*sigmaSqrd+tauSqrd*row$value)/varX
    posterior_var <- sigmaSqrd*tauSqrd/varX
    
    #hypothesis testing 
    bayesianHypothesisTesting(prior_mean, prior_var, posterior_mean,posterior_var, alpha)
  }
  return(cbind(dataTable, pHs))
}

rbind.with.progress <- function(progress, steps){  
  bayescount <- 0
  function(...) {
    new.entries <- ((length(list(...)) - 1) / steps) 
    bayescount <<- bayescount + new.entries
    if(!is.null(progress))
      progress(detail="Computing Bayes statistics", value=bayescount)            
    rbind(...)
  }
}

bayesianHypothesisTesting <- function(prior_mean, prior_var, posterior_mean, posterior_var, alpha=0.05){
  
  #It's the 0.05 and the 0.95 quantiles of the prior distribution
  #that will serve us as a threshold
  #for the posterior distribution
  prior_upper_quantile <- qnorm(1-alpha, mean=prior_mean, sd=sqrt(prior_var))
  prior_lower_quantile <- qnorm(alpha, mean=prior_mean, sd=sqrt(prior_var))
  
  #probability of alternative hypothesis H1: activation effect or P(?i - Theta0 > a|Xi)
  pH1 <- pnorm(prior_upper_quantile, mean=posterior_mean, sd=sqrt(posterior_var),lower.tail=F)
  
  #probability of alternative hypothesis H2: inhibition effect or P(?i - Theta0 < -a|Xi)
  pH2 <- pnorm(prior_lower_quantile, mean=posterior_mean, sd=sqrt(posterior_var), lower.tail=T)
  
  #probability of null hypothesis H0: no effect or P(|?i - Theta0| <= a|Xi)
  pH0 <- 1-(pH1+pH2)
  
  #probability of observing an effect
  pH3 <- pH1+pH2
  
  return(data.frame(p_no_effect=pH0, p_effect=pH3, p_promotor=pH1, p_suppressor=pH2))
}