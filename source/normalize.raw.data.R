normalizeRawData <- function(plates, control.based=F, pos.ctrl=NULL, neg.ctrl=NULL, updateProgress=NULL, compute.B=TRUE){
  library(dplyr)
  
  if(is.function(updateProgress)) updateProgress(detail="Plate based", value=0.2)
  #applying normalization strategies  
  plates.norm <- plates %>% group_by(Experiment, Readout, Plate, Replicate) %>% dplyr::mutate(
                       rzscore=(Raw - median(Raw, na.rm=T))/mad(Raw, na.rm=T),
                       zscore=(Raw - mean(Raw, na.rm=T))/sd(Raw, na.rm=T), 
                       centered=Raw/mean(Raw, na.rm=T),
                       rcentered=Raw/median(Raw, na.rm=T))
  
  if(control.based)
  {
    if(is.function(updateProgress)) updateProgress(detail="Control based", value=0.3)
       
    if(!is.null(neg.ctrl))
    {
      plates.norm <- plates.norm %>% dplyr::mutate(poc=(Raw / mean(Raw[Control==neg.ctrl], na.rm=T)))
    }
    if(!is.null(neg.ctrl) && !is.null(pos.ctrl)){
      plates.norm <- plates.norm %>% dplyr::mutate(npi=(mean(Raw[Control==neg.ctrl], na.rm=T) - Raw) 
                       / (mean(Raw[Control==neg.ctrl], na.rm=T) - mean(Raw[Control==pos.ctrl], na.rm=T)))
    }
  }                
  
  if(compute.B){
    if(is.function(updateProgress)) updateProgress(detail="Bscore", value=0.6)
    plates.norm <- do(plates.norm, Bscore(.))
  }
  #sort
  if(is.function(updateProgress)) updateProgress(detail="Sorting", value=0.8)
  plates.norm <- plates.norm %>% arrange(Experiment, Plate, Row, Column, Replicate) 
  
  plates.norm$wellCount <- as.integer(row.names(plates.norm))
  plates.norm$Replicate <- as.factor(plates.norm$Replicate)
  plates.norm$Plate <- as.factor(plates.norm$Plate)
  return(plates.norm)
}