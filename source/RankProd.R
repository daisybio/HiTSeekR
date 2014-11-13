require(stringr)
require(RankProd)
require(plyr)
require(reshape2)
require(ggplot2)

MCF7.matrix <- acast(MCF7, Plate + Well.position + Sample ~ Replicate, value.var="zscore")
MCF12.matrix <- acast(MCF12A, Plate + Well.position + Sample ~ Replicate, value.var="zscore")

rank.prod.matrix <- cbind(MCF7.matrix, MCF12.matrix)
colnames(rank.prod.matrix) <- c(rep("MCF7", 3), rep("MCF12", 3))
rank.prod.names <- row.names(MCF7.matrix)
#rank.prod.names <- sapply(str_split(row.names(MCF7.matrix), "_"), function(x){return(x[3])})

RP.out <- RP(rank.prod.matrix, c(rep(0,3), rep(1,3)), logged=F)

#killing MCF7, not-killing MCF12
result <- topGene(RP.out, cutoff=0.05, method="pval", logged=F, gene.names=rank.prod.names)[[1]]
result.matrix <- rank.prod.matrix[result[,1],]
MCF12A.mean <- apply(result.matrix[,c(1,2,3)], 1, mean)
result.matrix <- result.matrix[MCF12A.mean > quantile(MCF12A[,"zscore"], 0.1), MCF12A.mean < quantile(MCF12A[,"zscore"])]
plot.data <- melt(result.matrix)
colnames(plot.data) <- c("Sample", "Screen", "value")
plot.data$Sample <- reorder(plot.data$Sample, seq(1:nrow(plot.data)))
qplot(x=Sample, y=value, data=plot.data, fill=Screen, geom="boxplot", position="dodge") +  theme(axis.text.x=element_text(angle=45, hjust=1))