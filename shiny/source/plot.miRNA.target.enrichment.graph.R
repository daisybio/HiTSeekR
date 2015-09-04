library(foreach)
library(reshape2)
library(igraph)

prepare.kpm.output.for.plotting <- function(kpm.result, indicator.matrix, hit.list, screenType)
{
  #kpm.res <<- kpm.result
  #ind.m <<- indicator.matrix
  #hit.l <<- hit.list
  #browser()
  #get list of edges
  edges <- lapply(kpm.result$resultGraphs, function(x){return(lapply(x$edges, function(y){c(as.numeric(y$source), as.numeric(y$target))}))})
  
  if(length(edges) == 0) stop("No results were returned. Try different K and L settings.")
  #get union graph id
  unionGraph <- foreach(i = 1:length(kpm.result$resultGraphs)) %do%
  {
    if(kpm.result$resultGraphs[[i]]$isUnionSet) return(TRUE)
    else return(FALSE)
  }
  #extract edges
  unionGraph <- which(unlist(unionGraph))
  edges.df <- as.data.frame(t(as.data.frame(edges[[unionGraph]])))
  
  if(nrow(edges.df) == 0) stop("No results were returned. Try different K and L settings.")
  
  colnames(edges.df) <- c("source", "target")
  
  #remove edges that point to the node itself
  edges.df <- edges.df[-which(edges.df$source == edges.df$target),]
  
  #check if we have any edges left
  if(nrow(edges.df) == 0) stop("no edges found where source != target")
  
  nodes <- kpm.result$resultGraphs[[unionGraph]]$nodes
  
  #get node ids and overlap count
  node.ids <- foreach(x = nodes, .combine=rbind) %do%
  {
    data.frame('id'=x$id, 'overlapCount'=x$overlapCount) 
  }
  node.ids[,1] <- as.character(node.ids[,1])
  row.names(node.ids) <- node.ids[,1]
  
  if(screenType == "miRNA")
  {
    #get miRNA interactions
    target.interactions.matrix <- indicator.matrix[node.ids$id,]
    target.interactions <- target.interactions.matrix
    target.interactions$gene <- row.names(target.interactions)
    target.interactions <- melt(target.interactions)
    target.interactions <- subset(target.interactions, value==1)[,c(2,1)]
    colnames(target.interactions) <- c("source", "target")
    
    #add miRNA interactions to edge and node data frame
    edges.df <- rbind(edges.df, target.interactions)
    node.ids <- rbind(node.ids, data.frame(id=target.interactions$source, overlapCount=1))
  }
  return(list(edges.df, node.ids))
}

plot.miRNA.target.enrichment.graph.d3 <- function(kpm.result, indicator.matrix, hit.list, highlight=1, screenType)
{
  if(is.null(kpm.result))
  {
    return(NULL)
  }
  kpm.data <- prepare.kpm.output.for.plotting(kpm.result, indicator.matrix, hit.list, screenType)
  edges.df.mirna <- kpm.data[[1]]
  
  #default edge weight
  edges.df.mirna$value <- 1
  edges.df.levels <- unique(c(edges.df.mirna$source, edges.df.mirna$target))
  edges.df.mirna$source.id <- as.integer(factor(edges.df.mirna$source, levels=edges.df.levels))-1
  edges.df.mirna$target.id <- as.integer(factor(edges.df.mirna$target, levels=edges.df.levels))-1
  
  node.ids <- data.frame(id=factor(edges.df.levels, levels=edges.df.levels))
  
  #replace node gene ids with gene symbols
  #library(org.Hs.eg.db)
  #node.ids$symbol <- left_join(node.ids, as.data.frame(org.Hs.egSYMBOL), by=c("id" = "gene_id"))$symbol  
  
  if(screenType == "miRNA")
  {
    #color miRNAs according to suppressors / promoters
    categories <- left_join(node.ids, hit.list, by=c("id"="mature_name")) %>% group_by(id) %>% summarise(category=unique(category))
    node.ids <- left_join(node.ids, categories,by="id")
    node.ids[is.na(node.ids$category),"category"] <- "NA"
  }
  #node.ids$category <- as.integer(factor(node.ids$category, levels=c("promotor", "included", "gene", "suppressor")))
  
  #color genes found multiple times
  overlap <- left_join(node.ids, kpm.data[[2]], by="id") 
  
  node.ids[which(node.ids$id %in% unique(overlap[which(overlap$overlapCount > highlight), "id"])), "category"] <- "multiple"
  #add gene symbols
  symbols <- left_join(node.ids, as.data.frame(org.Hs.egSYMBOL), by=c("id" = "gene_id"))$symbol
  node.ids[!is.na(symbols), "id"] <- symbols[!is.na(symbols)]
  forceNetwork(Links = edges.df.mirna, Nodes = node.ids, Source="source.id", 
               Target="target.id", Value="value", NodeID = "id", Group="category",
               colourScale="d3.scale.ordinal().range(['#1f77b4', '#ff7f0e', '#2ca02c','#9b9e9b', '#d62728']).domain(['promotor', 'included', 'multiple', 'NA', 'suppressor']);")
}

plot.miRNA.target.enrichment.graph.igraph <- function(kpm.result, indicator.matrix, hit.list, screenType){
  
  if(is.null(kpm.result))
  {
    stop("Waiting for KPM results...")
  }
  kpm.data <- prepare.kpm.output.for.plotting(kpm.result, indicator.matrix, hit.list, screenType)
  edges.df.mirna <- kpm.data[[1]]
  node.ids <- kpm.data[[2]]

  #create graph
  g <- graph.data.frame(edges.df.mirna, directed=F)
  
  #add color gradient
  colGradientFunction <- colorRampPalette(c("grey", "green"))
  colGradient <- colGradientFunction(max(node.ids$overlapCount))
  node.ids$color <- colGradient[node.ids$overlapCount]
  
  #add colors (based on overlap count)
  g.colors <- node.ids[V(g)$name,"color"]
  g.colors[is.na(g.colors)] <- "#87898a"
  V(g)$color <- g.colors
  
  #return(list(g, target.interactions))
  
  #replace node gene ids with gene symbols
  library(org.Hs.eg.db)
  node.names <- V(g)$name
  symbols <- left_join(data.frame(name = node.names), as.data.frame(org.Hs.egSYMBOL), by=c("name" = "gene_id"))$symbol
  V(g)$name[-which(is.na(symbols))] <- symbols[-which(is.na(symbols))]
  
  #color miRNAs according to suppressors / promoters
  if(screenType == "miRNA"){
    category <- left_join(data.frame(name = node.names), hit.list, by=c("name"="mature_name"))$category
    V(g)$color[which(category == "promotor")] <- "#80B1D3"
    V(g)$color[which(category == "suppressor")] <- "#FB8072"
    V(g)$color[which(category == "included")] <- "#FDB462"
  }
  else if(screenType == "siRNA"){
    category <- left_join(data.frame(name = node.names), hit.list, by=c("name" = "gene_id"))$category
    V(g)$color[which(category == "promotor")] <- "#80B1D3"
    V(g)$color[which(category == "suppressor")] <- "#FB8072"
    V(g)$color[which(category == "included")] <- "#FDB462"
  }
  
  #layout
  layout <- layout.fruchterman.reingold(g)
  
  #plot graph
  plot(g, layout=layout)
  
  if(screenType == "miRNA"){
    legend("right", legend=c("Included", "Suppressor", "Promoter", "Reocurring gene"), fill=c("#FDB462", "#FB8072", "#80B1D3", "green"))
  }
}
