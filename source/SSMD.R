ssmd <- function(exp.data, neg.ctrl, signalColumn, summarise.results=TRUE, updateProgress)
{        
  #update progress bar
  updateProgress(detail="Computing SSMD scores", value=ssmdPlateCounter/ssmdPlateMax)
  ssmdPlateCounter <<- ssmdPlateCounter + 1
  
  neg.ctrls <- dplyr::filter(exp.data, Control %in% neg.ctrl)
  if(nrow(neg.ctrls) == 0) stop("At least one plate is missing negative controls")
  result <- foreach(nc = neg.ctrl, .combine=rbind) %do%
  {      
    neg.ctrl.data <- exp.data %>% dplyr::filter(Control == nc)
    if(nrow(neg.ctrl.data) < 3) stop("Cannot estimate variance in SSMD with < 3 negative control wells per plate")
    neg.ctrl.mean <- mean(neg.ctrl.data[[signalColumn]], na.rm=T)
    neg.ctrl.sd <- sd(neg.ctrl.data[[signalColumn]], na.rm=T)
    
    calc.ssmd <- function(x)
    {      
      sampleSignal <- x[[signalColumn]]
      
      if(nrow(x) == 1)
      {        
        return ((sampleSignal - neg.ctrl.mean) / (sqrt(2) * neg.ctrl.sd))
      }
      else
      {
        return ((mean(sampleSignal, na.rm=T) - neg.ctrl.mean) / sqrt((neg.ctrl.sd)^2 + sd(sampleSignal, na.rm=T)^2))    
      }
    }    
    result <- as.data.frame(exp.data)
    result <- result %>% dplyr::group_by(Sample) %>% 
      dplyr::do(NEG.CTRL=nc, SSMD=calc.ssmd(.))       
    return(na.omit(result))
  }
  if(summarise.results)
    result <- result %>% dplyr::group_by(Sample) %>% dplyr::summarise(SSMD = mean(unlist(SSMD), na.rm=T))
  else{
    result <- result %>% dplyr::mutate(SSMD = unlist(SSMD), NEG.CTRL = unlist(NEG.CTRL))
  }
  return(as.data.frame(result))
}