#heatmap
my.heatmap <- function(plates, signalType, margin, method, outliers, ncol=3, title="", showSampleLabels=T, signalColumn="Raw")
{
  kIn1.outliers <- outliers
  
  require(ggplot2)
  require(grid)
  #title <- paste(title, signalType ,": Heatmaps (labels for +-", margin, method, ")")
  p2 <- qplot(Column, Row, data=plates, xlab="Column", ylab="Row")
  p2 <- p2 + geom_raster(aes_string(fill = signalColumn));
  p2 <- p2 + scale_fill_gradient2(midpoint=mean(plates[,signalColumn], na.rm=T), mid="#FFFFCC")
  #p2 <- p2 + scale_fill_gradient(low = colorA, high = colorB, name=paste(signalType));
  p2 <- p2 + facet_wrap(~Plate, ncol=ncol);
  
  if(dim(kIn1.outliers)[1] > 0)
  {
    bottomRight <- subset(kIn1.outliers,Row>1 & Column > 1)
    
    if(dim(bottomRight)[1] > 0)
    {
      if(showSampleLabels){
        p2 <- p2 + geom_text(data=bottomRight, aes(label=Sample), vjust=-1.5, hjust=.7)
      }
      p2 <- p2 + geom_segment(data=bottomRight, aes(x=Column-0.5, y=Row-0.5, xend=Column, yend=Row), colour=I("black"), arrow=arrow(angle=45, length=unit(0.2, "cm")))
    }
    
    firstRow <- subset(kIn1.outliers, Row==1)
    
    if(dim(firstRow)[1] > 0)
    {
      evenColumns <- subset(firstRow, (Column %% 2) == 0)
      if(dim(evenColumns)[1] > 0 )
      {
        if(showSampleLabels){
          p2 <- p2 + geom_text(data=evenColumns, aes(label=Sample), vjust=2.0, hjust=0)
        }
        p2 <- p2 + geom_segment(data=evenColumns, aes(x=Column+0.5, y=Row+0.5, xend=Column, yend=Row), colour=I("black"), arrow=arrow(angle=45, length=unit(0.2, "cm")))
      }
      
      unevenColumns <- subset(firstRow, (Column %% 2) == 1)
      if(dim(unevenColumns)[1] > 0)
      {
        if(showSampleLabels){
          p2 <- p2 + geom_text(data=unevenColumns, aes(label=Sample), vjust=3.0, hjust=0)
        }
        p2 <- p2 + geom_segment(data=unevenColumns, aes(x=Column+1.0, y=Row+1.0, xend=Column, yend=Row), colour=I("black"), arrow=arrow(angle=45, length=unit(0.2, "cm")))
      }
    }
    
    firstColumn <- subset(kIn1.outliers, Row>1 & Column == 1)
    
    if(dim(firstColumn)[1] >0 )
    {
      if(showSampleLabels){
        p2 <- p2 + geom_text(data=firstColumn, aes(label=Sample), vjust=-1.5, hjust=0.3)
      }
      p2 <- p2 + geom_segment(data=firstColumn, aes(x=Column+0.5, y=Row-0.5, xend=Column, yend=Row), colour=I("black"), arrow=arrow(angle=45, length=unit(0.2, "cm")))
    }
  }
  
  p2 <- p2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.margin=unit(0.1, "lines"), panel.margin=unit(0, "lines"), plot.margin=unit(c(1, 1, 0.5, 0.5), "lines"), 
                  plot.title=element_text(size=18), strip.background=element_rect(fill="grey90", colour="grey50"))
  #p2 <- p2 + geom_dl(aes(label=Sample), method=list("smart.grid"), data=kIn1.outliers)
  p2 <- p2 + scale_x_continuous(expand=c(0,0), breaks=seq(2,length(unique(plates$Plate)),3));
  p2 <- p2 + scale_y_reverse(expand=c(0,0), breaks=seq(2,length(unique(plates$Plate)),3));
  print(p2);
}
