normalizeRawData <- function(plates, control.based=F, pos.ctrl=NULL, neg.ctrl=NULL){
  library(dplyr)
  #applying normalization strategies
  plates.norm <- plates %>% group_by(Experiment, Plate, Replicate) %>% mutate(
                       rzscore=(Raw - median(Raw, na.rm=T))/mad(Raw, na.rm=T),
                       zscore=(Raw - mean(Raw, na.rm=T))/sd(Raw, na.rm=T), 
                       centered=Raw/mean(Raw, na.rm=T),
                       rcentered=Raw/median(Raw, na.rm=T))
  if(control.based)
  {
    if(!is.null(neg.ctrl))
    {
      plates.norm <- plates.norm %>% mutate(poc=(Raw / mean(Raw[Control==neg.ctrl], na.rm=T)))
    }
    if(!is.null(neg.ctrl) && !is.null(pos.ctrl)){
      plates.norm <- plates.norm %>% mutate(npi=(mean(Raw[Control==neg.ctrl], na.rm=T) - Raw) 
                       / (mean(Raw[Control==neg.ctrl], na.rm=T) - mean(Raw[Control==pos.ctrl], na.rm=T)))
    }
  }                
  plates.norm <- do(plates.norm, Bscore(.))
  #library(plyr)
  #plates.norm <- ddply(plates.norm, .(Plate, Replicate), Bscore)
  #plates.norm <- posEffectNorm(plates.norm)
  #plate comparison
  #sort
  plates.norm <- plates.norm %>% arrange(Experiment, Plate, Row, Column, Replicate) 
  
  plates.norm$wellCount <- as.integer(row.names(plates.norm))
  plates.norm$Replicate <- as.factor(plates.norm$Replicate)
  plates.norm$Plate <- as.factor(plates.norm$Plate)
  return(plates.norm)
}