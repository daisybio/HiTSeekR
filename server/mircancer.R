# load database #
mircancer.database <- reactive({
  read.table("data/miRCancerJune2013.txt", sep="\t", header=T, quote="\"")
})

# find entries for currently selected hits #
outliers.mircancer <- reactive({
  outliers <- outliers()
  mirdb <- mircancer.database()
  outliers$Sample <- sub("hsa-miR", "hsa-mir", outliers$Sample)
  outliers <- outliers[-grep(".*>NA</a>", outliers$Accession),]
  merge(outliers[,c("Sample", "Accession")], mirdb, by.x="Sample", by.y="mirId", all.x=T)
})