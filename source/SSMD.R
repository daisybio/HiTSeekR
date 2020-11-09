#TODO: change the name of variables that are called ssmd to quality control
ssmd <- function(exp.data, neg.ctrl, signalColumn, summarise.results=TRUE, updateProgress, use.zfactor = F)
{        
  quality.control.string <- ifelse(use.zfactor, "Z-factor", "SSMD")
  #update progress bar
  updateProgress(detail=sprintf("Computing %s scores", quality.control.string), value=ssmdPlateCounter/ssmdPlateMax)
  ssmdPlateCounter <<- ssmdPlateCounter + 1
  
  neg.ctrls <- dplyr::filter(exp.data, Control %in% neg.ctrl)
  
  if(nrow(neg.ctrls) == 0) stop("At least one plate is missing negative controls")
  result <- foreach(nc = neg.ctrl, .combine=rbind) %do%
  {      
    neg.ctrl.data <- exp.data %>% dplyr::filter(Control == nc)
    #TODO calculate based on all sample wells?
    if(nrow(neg.ctrl.data) < 3) stop(sprintf("Cannot estimate variance in %s with < 3 negative control wells per plate", 
                                             quality.control.string))
    neg.ctrl.mean <- mean(neg.ctrl.data[[signalColumn]], na.rm=T)
    neg.ctrl.sd <- sd(neg.ctrl.data[[signalColumn]], na.rm=T)
    
    calc.ssmd <- function(x, zfactor.instead)
    {      
      sampleSignal <- x[[signalColumn]]
      if(!zfactor.instead){
        if(nrow(x) == 1)
        {        
          return ((sampleSignal - neg.ctrl.mean) / (sqrt(2) * neg.ctrl.sd))
        }
        else
        {
          return ((mean(sampleSignal, na.rm=T) - neg.ctrl.mean) / sqrt((neg.ctrl.sd)^2 + sd(sampleSignal, na.rm=T)^2))    
        }
      } else {
        return (1 - (3 * (sd(sampleSignal, na.rm = T) + neg.ctrl.sd) 
                     / (abs(mean(sampleSignal, na.rm = T) - neg.ctrl.mean))))
      }
    }    
    result <- as.data.frame(exp.data)
    result <- result %>% dplyr::group_by(Sample) %>% 
      dplyr::do(NEG.CTRL=nc, SSMD=calc.ssmd(., use.zfactor))       
    return(na.omit(result))
  }
  if(summarise.results)
    result <- result %>% dplyr::group_by(Sample) %>% dplyr::summarise(SSMD = mean(unlist(SSMD), na.rm=T))
  else{
    result <- result %>% dplyr::mutate(SSMD = unlist(SSMD), NEG.CTRL = unlist(NEG.CTRL))
  }
  return(as.data.frame(result))
}