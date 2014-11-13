# load database #
mircancer.database <- reactive({
  read.table("data/miRCancerJune2013.txt", sep="\t", header=T, quote="\"")
})

# find entries for currently selected hits #
outliers.mircancer <- reactive({
  outliers <- outliers()
  mirdb <- mircancer.database()
  outliers$miRBase.ID.miRPlus.ID <- sub("hsa-miR", "hsa-mir", outliers$miRBase.ID.miRPlus.ID)
  outliers <- outliers[-grep(".*>NA</a>", outliers$miRBase.accession),]
  merge(outliers[,c("miRBase.ID.miRPlus.ID", "miRBase.accession")], mirdb, by.x="miRBase.ID.miRPlus.ID", by.y="mirId", all.x=T)
})