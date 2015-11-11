#In this function the actual hit detection is done.
find.hits.call <- function(exp.data, rep.data, method, margin, neg.ctrl, signalColumn, updateProgress, upperCutoff, lowerCutoff){
  #percentage cutoff is the simples method
  if(method == "cutoff")
  {
    outl <- dplyr::filter_(exp.data, interp("signalColumn > uC | signalColumn < lC", .values=list(signalColumn = as.name(signalColumn), uC=upperCutoff, lC=lowerCutoff)))
  }
  #Bayesian hit detection using negative controls and plate means to calculate priors.
  #Sample effect is then calculated as posterior probability
  else if(method == "Bayes")
  {    
    outl <- bayesianHitSelection(exp.data, neg.ctrl=neg.ctrl, signalColumn=signalColumn,alpha = 0.05, updateProgress=updateProgress)
  } 
  else if(margin == 0) outl <- exp.data #if we don't use a detection margin input == output
  else if(method == "SSMD") #SSMD needs to be treated separately because its based on replicates while the other methods operate on merged replicates
  {        
    result <- rep.data %>% dplyr::group_by(Plate) 
    ssmdPlateCounter <<- 0
    ssmdPlateMax <<- length(dplyr::group_size(result))
    result <- result %>% do(ssmd(., neg.ctrl, signalColumn, updateProgress=updateProgress))
    outl <- exp.data
    outl <- dplyr::left_join(outl,result, by=c("Plate", "Sample"))
    outl <- outl %>% dplyr::filter(abs(SSMD) >= margin)    
  }
  else{ #for all other methods, e.g. SD, MD and quartile
    outl <- find.hits(exp.data, method, margin, signalColumn=signalColumn, updateProgress=updateProgress)    
  }  
  
  outl <- as.data.frame(outl)
  #no hits were found? bad parameters
  if(nrow(outl) == 0) stop("No hits were found with these settings.")
  
  #add category labels to the samples. We naively use the mean of all signal values as middle point.
  outl[which(outl[,signalColumn] > mean(exp.data[,signalColumn], na.rm=T)),"category"] <- "promotor"
  outl[which(outl[,signalColumn] < mean(exp.data[,signalColumn], na.rm=T)),"category"] <- "suppressor"
  
  #return results
  return(outl)
}
