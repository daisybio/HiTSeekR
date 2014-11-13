processJos <- function(jos_raw)
{
  library(foreach)
  
  firstRow = 1
  columnSkip = 17
  currentCol = 1
  
  foreach(plate=1:18, .combine=rbind) %do%
  {
    foreach(replicate=1:2, .combine=rbind) %do% 
    {
      result <- foreach(column=1:12, .combine=rbind) %do%
      {
        foreach(row=1:8, .combine=rbind) %do%
        {
          #debug output
          #print(jos_raw[firstRow + row, currentCol + column])
          data.frame(Plate=plate, Row=row, Column=column, Replicate=replicate, Raw=jos_raw[firstRow + row, currentCol + column])    
        }
      }
      currentCol = currentCol + columnSkip 
      return(result)
    }
  }
}

josB <- function(plates){
  plates <- ddply(plates, .(Replicate, Plate, Row), function(x){ x$rowMedian <-median(x$Raw, na.rm=T); return(x) })
  plates <- ddply(plates, .(Plate, Replicate), 
                  function(x){ x$josB <- 
                                 ((x$Raw - median(x$Raw * median(x$Raw, na.rm=T) / x$rowMedian, na.rm=T))
                                  /mad(x$Raw, na.rm=T, constant=1)); 
                               return(x)})
  plates <- ddply(plates, .(), transform, josB = round(josB, 3))
  return(plates[,-(ncol(plates)-1)])
}
jos_raw <- read.delim("jos_raw.txt")
jos_processed <- processJos(jos_raw)

#add Jos Bscore (only row normalization)
jos_processed <- josB(jos_processed)

#read sample names with Bscores from paper
jos_names <- read.delim("A375_MTS_names.txt")
jos_names <- ddply(jos_names, .(), transform, MTS.1=round(MTS.1,3), MTS.2=round(MTS.2,3))

#merge our data frame with sample names from plos one paper suppl. file
jos_mirnas <- merge(jos_processed, jos_names[,c(2,3)], by.x="josB", by.y="MTS.1")[,c(3,4,5,6,8)]
jos_mirnas <- rbind(jos_mirnas, merge(jos_processed, jos_names[,c(2,4)], by.x="josB", by.y="MTS.2")[,c(3,4,5,6,8)])
colnames(jos_mirnas)[5] <- "Sample"
jos_processed <- merge(jos_processed, jos_mirnas, by=c("Plate", "Row", "Column", "Replicate"), all.x=T)
write.table(jos_processed, "jos_processed.txt", row.names=F, sep="\t", quote=F)
