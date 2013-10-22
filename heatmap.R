#outlier detection
my.outliers <- function(plates, method, margin=2, withControls=F)
{ 
  upper_margin <- margin
  
  if(withControls)
  {
    data <- plates
  }
  else
  {
    data <- subset(plates, column < 12)
  }
  
  if(method=="SD")
  {
    kIn1.mean <- mean(plates$signal, na.rm=T)
    kIn1.sd <- sd(plates$signal, na.rm=T)
    kIn1.outliers <-  subset(data,((signal > (kIn1.mean + (upper_margin * kIn1.sd))) 
                                   | (signal < (kIn1.mean - (margin * kIn1.sd))))) 
  }
  else if(method=="MAD")
  {
    kIn1.median <- median(plates$signal, na.rm=T)
    kIn1.mad <- mad(plates$signal, na.rm=T)
    kIn1.outliers <- subset(data,((signal > (kIn1.median + (upper_margin * kIn1.mad)))
                                  | (signal < (kIn1.median - (margin * kIn1.mad)))))
  }
  
  else if(method=="quartile")
  {
    kIn1.quantiles <- quantile(plates$signal, na.rm=T)
    kIn1.IQR <- kIn1.quantiles[4] - kIn1.quantiles[2]
    kIn1.outliers <- subset(data, ((signal > (kIn1.quantiles[4] + (upper_margin * kIn1.IQR))) 
                                   | (signal < (kIn1.quantiles[2] - (margin * kIn1.IQR)))))
  }

  else kIn1.outliers <- plates
  
  return(kIn1.outliers)
}

#heatmap
my.heatmap <- function(plates, signalType, margin, method, outliers, ncol=3, title="", colorA="red", colorB="yellow")
{
  kIn1.outliers <- outliers
  
  require(ggplot2)
  require(grid)
  title <- paste(title, signalType ,": Heatmaps (labels for +-", margin, method, ")")
  p2 <- qplot(column, row, data=plates, main=title, xlab="column", ylab="row")
  p2 <- p2 + geom_tile(line=0, aes(fill = signal));
  #p2 <- p2 + scale_fill_gradient2(midpoint=mean(kIn$log10, na.rm=T));
  p2 <- p2 + scale_fill_gradient(low = colorA, high = colorB, name=paste(signalType));
  p2 <- p2 + facet_wrap(~Plate, ncol=ncol);
  
  if(dim(kIn1.outliers)[1] > 0)
  {
    bottomRight <- subset(kIn1.outliers,row>1 & column > 1)
    
    if(dim(bottomRight)[1] > 0)
    {
      p2 <- p2 + geom_text(data=bottomRight, aes(label=Sample), vjust=-1.5, hjust=.7)
      p2 <- p2 + geom_segment(data=bottomRight, aes(x=column-0.5, y=row-0.5, xend=column, yend=row), colour=I("black"), arrow=arrow(angle=45, length=unit(0.2, "cm")))
    }
    
    firstRow <- subset(kIn1.outliers, row==1)
    
    if(dim(firstRow)[1] > 0)
    {
      evenColumns <- subset(firstRow, (column %% 2) == 0)
      if(dim(evenColumns)[1] > 0 )
      {
        p2 <- p2 + geom_text(data=evenColumns, aes(label=Sample), vjust=2.0, hjust=0)
        p2 <- p2 + geom_segment(data=evenColumns, aes(x=column+0.5, y=row+0.5, xend=column, yend=row), colour=I("black"), arrow=arrow(angle=45, length=unit(0.2, "cm")))
      }
      
      unevenColumns <- subset(firstRow, (column %% 2) == 1)
      if(dim(unevenColumns)[1] > 0)
      {
        p2 <- p2 + geom_text(data=unevenColumns, aes(label=Sample), vjust=3.0, hjust=0)
        p2 <- p2 + geom_segment(data=unevenColumns, aes(x=column+1.0, y=row+1.0, xend=column, yend=row), colour=I("black"), arrow=arrow(angle=45, length=unit(0.2, "cm")))
      }
    }
    
    firstColumn <- subset(kIn1.outliers, row>1 & column == 1)
    
    if(dim(firstColumn)[1] >0 )
    {
      p2 <- p2 + geom_text(data=firstColumn, aes(label=Sample), vjust=-1.5, hjust=0.3)
      p2 <- p2 + geom_segment(data=firstColumn, aes(x=column+0.5, y=row-0.5, xend=column, yend=row), colour=I("black"), arrow=arrow(angle=45, length=unit(0.2, "cm")))
    }
  }
  
  p2 <- p2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.margin=unit(0.1, "lines"), panel.margin=unit(0, "lines"), plot.margin=unit(c(1, 1, 0.5, 0.5), "lines"), 
                  plot.title=element_text(size=18), strip.background=element_rect(fill="grey90", colour="grey50"))
  #p2 <- p2 + geom_dl(aes(label=Sample), method=list("smart.grid"), data=kIn1.outliers)
  p2 <- p2 + scale_x_continuous(expand=c(0,0), breaks=seq(2,12,3));
  p2 <- p2 + scale_y_reverse(expand=c(0,0), breaks=seq(2,12,3));
  print(p2);
}
