#hit detection
find.hits <- function(plates, method, margin=2, withControls=F, signalColumn="Raw", updateProgress=NULL)
{ 
  library(dplyr)
  
  if(!is.null(updateProgress)) updateProgress(detail=method)
  
  upper_margin <- margin
  
  if(withControls)
  {
    data <- plates
  }
  else
  {
    data <- filter(plates, is.na(Control))
  }
  if(method=="SD")
  {
    kIn1.mean <- mean(data[,signalColumn], na.rm=T)
    kIn1.sd <- sd(data[,signalColumn], na.rm=T)
    upper_limit <- kIn1.mean + (upper_margin * kIn1.sd)
    lower_limit <- kIn1.mean - (margin * kIn1.sd)
  }
  else if(method=="MAD")
  {
    kIn1.median <- median(data[,signalColumn], na.rm=T)
    kIn1.mad <- mad(data[,signalColumn], na.rm=T)
    upper_limit <- kIn1.median + (upper_margin * kIn1.mad)
    lower_limit <- kIn1.median - (margin * kIn1.mad)
  }
  
  else if(method=="quartile")
  {
    kIn1.quantiles <- quantile(data[,signalColumn], na.rm=T)
    kIn1.IQR <- kIn1.quantiles[4] - kIn1.quantiles[2]
    upper_limit <- kIn1.quantiles[3] + (upper_margin * kIn1.IQR)
    lower_limit <- kIn1.quantiles[3] - (margin * kIn1.IQR)    
  }
  
  else return(plates)
  
  result <- data %>% filter(get(signalColumn, envir=as.environment(data)) > upper_limit 
                            | get(signalColumn, envir=as.environment(data)) < lower_limit)
  return(result)
}
