library(foreach)
library(reshape2)
library(igraph)
library(org.Hs.eg.db)

prepare.kpm.output.for.plotting <- function(kpm.result, indicator.matrix, hit.list, screenType, selectedGraph)
{
  #kpm.res <<- kpm.result
  #ind.m <<- indicator.matrix
  #hit.l <<- hit.list
  
  #get list of edges
  edges <- lapply(kpm.result$resultGraphs, function(x){return(lapply(x$edges, function(y){c(as.numeric(y$source), as.numeric(y$target))}))})
  
  if(length(edges) == 0) stop("No results were returned. Try different K and L settings.")
  #get union graph id
  result.selectedGraph <- foreach(i = 1:length(kpm.result$resultGraphs)) %do%
  {
    if(kpm.result$resultGraphs[[i]]$isUnionSet && is.na(selectedGraph)) return(TRUE)
    else if(!is.na(selectedGraph) && kpm.result$resultGraphs[[i]]$nodeSetNr == selectedGraph) return(TRUE)
    else return(FALSE)
  }
  #extract edges
  result.selectedGraph <- which(unlist(result.selectedGraph))
  edges.df <- as.data.frame(t(as.data.frame(edges[[result.selectedGraph]])))
  
  if(nrow(edges.df) == 0) stop("No results were returned. Try different K and L settings.")
  
  colnames(edges.df) <- c("source", "target")
  
  #remove edges that point to the node itself
  self_loops <- which(edges.df$source == edges.df$target)
  if(length(self_loops) > 0)
    edges.df <- edges.df[-which(edges.df$source == edges.df$target),]
  
  nodes <- kpm.result$resultGraphs[[result.selectedGraph]]$nodes
  
  #get node ids and overlap count
  node.ids <- foreach(x = nodes, .combine=rbind) %do%
  {
    data.frame('id'=x$id, 'overlapCount'=x$overlapCount) 
  }
  node.ids[,1] <- as.character(node.ids[,1])
  row.names(node.ids) <- node.ids[,1]
  
  if(screenType %in% c("miRNA", "compound"))
  {
    #get miRNA interactions
    target.interactions.matrix <- indicator.matrix[node.ids$id,]
    target.interactions <- target.interactions.matrix
    target.interactions$gene <- row.names(target.interactions)
    target.interactions <- melt(target.interactions)
    target.interactions <- subset(target.interactions, value==1)[,c(2,1)]
    colnames(target.interactions) <- c("source", "target")
    
    #add interactions to edge and node data frame
    edges.df <- rbind(edges.df, target.interactions)
    node.ids <- rbind(node.ids, data.frame(id=target.interactions$source, overlapCount=1))
  }
  
  return(list(edges.df, node.ids))
}

prepare.kpm.dynamic <- function(kpm.data, hit.list, screenType)
{
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
  }
  else if(screenType == "siRNA")
  {
    #color siRNAs according to suppressors / promoters
    categories <- left_join(node.ids, hit.list, by=c("id"="gene_id")) %>% group_by(id) %>% summarise(category=unique(category))
    node.ids <- left_join(node.ids, categories,by="id")
    node.ids[is.na(node.ids$category),"category"] <- ""
  }
  else if(screenType == "compound"){
    #color compounds according to suppressors / promoters
    categories <- left_join(node.ids, hit.list, by=c("id"="PubChem_CID")) %>% group_by(id) %>% summarise(category=unique(category))
    node.ids <- left_join(node.ids, categories,by="id")
  }
  
  #color genes found multiple times
  overlap <- kpm.data[[2]] %>% group_by(id) %>% summarise(overlapCount = max(overlapCount))
  node.ids <- left_join(node.ids, overlap, by="id") 
  
  #add gene symbols
  symbols <- left_join(node.ids, as.data.frame(org.Hs.egSYMBOL), by=c("id" = "gene_id"))$symbol
  node.ids[!is.na(symbols), "id"] <- symbols[!is.na(symbols)]
  node.ids[which(node.ids$category == ""), "category"] <- "gene"
  return(list(nodes = node.ids, edges = edges.df.mirna))
}


plot.kpm.drugst.one <- function(kpm.data, hit.list, screenType)
{
  plot_data <- prepare.kpm.dynamic(kpm.data, hit.list, screenType)
  
  edges <- data.frame(from = plot_data$nodes[plot_data$edges$source.id + 1, 'id'],
                      to = plot_data$nodes[plot_data$edges$target.id + 1, 'id'])
  plot_data$edges <- edges
  plot_data$nodes <- S4Vectors::rename(plot_data$nodes, category = 'group')
  plot_data$nodes$overlapCount <- NULL
  nodeGroups <- list(promotor = list(groupName = 'promotor', color = '#1f77b4', shape = 'circle', type = 'promotor'),
                     included = list(groupName = 'included', color = '#ff7f0e', shape = 'circle', type = 'included'),
                     gene = list(groupName = 'gene', color = '#9b9e9b', shape = 'circle', type = 'gene'),
                     suppressor = list(groupName = 'suppressor', color = '#d62728', shape = 'circle', type = 'suppressor'))
  config <- list(identifier = "symbol", nodeShadow = FALSE, edgeShadow = FALSE, physicsOn = TRUE, nodeGroups = nodeGroups)
  HTML(sprintf('<drugst-one
          id=\'drugstone-component-id\'
          config=\'%s\'
          network=\'%s\'>
        </drugst-one>', 
               jsonlite::toJSON(config, auto_unbox = TRUE),
               jsonlite::toJSON(plot_data)))
}


plot.kpm.d3 <- function(kpm.data, hit.list, screenType)
{
  plot_data <- prepare.kpm.dynamic(kpm.data, hit.list, screenType)
  
  forceNetwork(Links = plot_data$edges, Nodes = plot_data$nodes, Nodesize="overlapCount", 
               Source="source.id", opacity=0.9,legend = TRUE, opacityNoHover=0.8,
               Target="target.id", Value="value", NodeID = "id", Group="category", zoom = TRUE,
               colourScale="d3.scale.ordinal().range(['#1f77b4', '#ff7f0e', '#9b9e9b', '#d62728']).domain(['promotor', 'included', 'gene', 'suppressor']);")
}

plot.kpm.igraph <- function(kpm.data, hit.list, screenType){
  
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
  node.names <- V(g)$name
  symbols <- left_join(data.frame(name = node.names, stringsAsFactors = FALSE), as.data.frame(org.Hs.egSYMBOL), by=c("name" = "gene_id"))$symbol
  if(any(is.na(symbols)))
  {
    V(g)$name[-which(is.na(symbols))] <- symbols[-which(is.na(symbols))]
  }
  else{
    V(g)$name <- symbols
  }
  
  #color miRNAs according to suppressors / promoters
  if(screenType == "miRNA"){
    category <- left_join(data.frame(name = node.names, stringsAsFactors = FALSE), hit.list, by=c("name"="mature_name"))$category
    V(g)$color[which(category == "promotor")] <- "#80B1D3"
    V(g)$color[which(category == "suppressor")] <- "#FB8072"
    V(g)$color[which(category == "included")] <- "#FDB462"
  }
  else if(screenType == "siRNA"){
    category <- left_join(data.frame(name = node.names, stringsAsFactors = FALSE), hit.list, by=c("name" = "gene_id"))$category
    V(g)$color[which(category == "promotor")] <- "#80B1D3"
    V(g)$color[which(category == "suppressor")] <- "#FB8072"
    V(g)$color[which(category == "included")] <- "#FDB462"
  }
  else if(screenType == "compound"){
    category <- left_join(data.frame(name = node.names, stringsAsFactors = FALSE), hit.list, by=c("name" = "PubChem_CID"))$category
    V(g)$color[which(category == "promotor")] <- "#80B1D3"
    V(g)$color[which(category == "suppressor")] <- "#FB8072"
    V(g)$color[which(category == "included")] <- "#FDB462"
  }
  
  #layout
  layout <- layout.fruchterman.reingold(g)
  
  #plot graph
  plot(g, layout=layout, margin=c(0,-0.5,0,0))
  
  if(screenType == "miRNA"){
    legend("right", legend=c("Included", "Suppressor", "Promoter", "Reocurring gene"), fill=c("#FDB462", "#FB8072", "#80B1D3", "green"))
  }
  else if(screenType == "siRNA"){
    legend("right", legend=c("Included", "Suppressor", "Promoter"), fill=c("#FDB462", "#FB8072", "#80B1D3"))
  }
}
