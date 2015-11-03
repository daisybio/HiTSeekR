posEffectNorm <- function(plates, Rows=T, Cols=T)
{
  plates.norm <- ddply(plates, .(Plate, Replicate), transform, rcentered=Raw/median(Raw, na.rm=T))
  
  if(Rows){
    rE <- ddply(plates.norm, .(Row), summarise, rE = median(rcentered, na.rm=T))$rE
    #rE <- 1 + (rE - median(rE))
  } 
  else rE <- ddply(plates.norm, .(Row), summarise, rE = 1)$rE
  
  if(Cols){
    cE <- ddply(plates.norm, .(Column), summarise, cE = median(rcentered, na.rm=T))$cE
    #cE <- 1 + (cE - median(cE))
  }
  else cE <- ddply(plates.norm, .(Column), summarise, cE = 1)$cE
     
  calcPosEffectNorm <- function(x, rE, cE)
  {
    return(x$Raw / (rE[x$Row] * cE[x$Column]))
  }
  
  result <- adply(plates, 1, calcPosEffectNorm, rE=rE, cE=cE)                
  colnames(result)[length(colnames(result))] <- "posEffectNorm"
  result <- ddply(result, .(Plate, Replicate), transform, posEffectNorm=((posEffectNorm-median(posEffectNorm, na.rm=T))/mad(posEffectNorm, na.rm=T)))
  return(result)
}