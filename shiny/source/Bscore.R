Bscore <-function(plate, value.var="Raw"){
  
  require(reshape2)
  require(plyr)
  #casting the signal to a matrix
  values <- acast(plate, Row~Column, value.var=value.var)
  #Tukey's two way median polish
  med.pol <- medpolish(values, trace.iter=F, na.rm=T)
  
  #extract values
  estPlateAverage <- med.pol$overall
  rowEffects <- med.pol$row
  colEffects <- med.pol$col
  
  #computation function
  calcB <- function(x, mu_p, rE, cE){ return(x[[value.var]] - (mu_p + rE[x$Row] + cE[x$Column]))}
  #calculate Bscore for each well
  result <- adply(plate,1, calcB, mu_p=estPlateAverage, rE=rowEffects, cE=colEffects)
  colnames(result)[ncol(result)] <- "Bscore"      
  result$Bscore <- result$Bscore / mad(result[[value.var]], na.rm=T)
  return(result)
}