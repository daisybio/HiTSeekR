## find mRNA targets ##
targetsForInteractionGraph <- reactive({
  
  input$updateTargets
  
  if(input$useConsensus == "hit list") data <- outliers()
  else data <- consensusHitList()
  
  result <- isolate({    
    result <- data.frame()
    grouped.targets <- getTargets(data, group.miRNAs=input$group.miRNAs, group.miRNAs.threshold=input$group.miRNAs.threshold, databases=input$selectedTargetDBs)  

    if(input$remove.nas.from.target.list){
      grouped.targets <- na.omit(grouped.targets)
    }
    
    grouped.targets$gene_symbol <- paste("<a href='http://www.ncbi.nlm.nih.gov/gene/?term=", 
                                       grouped.targets$gene_id, "%5Buid%5D", "' target='_blank'>",
                                       grouped.targets$gene_symbol,
                                       "</a>", sep="")
    

#     for(i in 1:nrow(grouped.targets)){
#       temp.result <- data.frame(
#         gene_symbol <- gsub("<|>", "", str_extract(grouped.targets[i, "gene_symbol"], ">.*?<")), 
#         miRNA <- str_split(grouped.targets[i, "miRNA_list"], "/")
#       )
#       
#       colnames(temp.result) <- c("gene_symbol", "miRNA")
#       result <- rbind(result, temp.result)
#     }
#     data$Sample <- sub("mir", "miR", data$Sample)
#     result <- merge(result, unique(data[,c("Sample", "category")]), by.x="miRNA", by.y="Sample", all.x=T, all.y=F)
#     result[,1] <- gsub("hsa-", "", result[,1])
    #print(result)
    return(grouped.targets)
  })
  
  return(result)
})

# formatted targets #
mirna.targets <- reactive({
  
  input$updateTargets
  
  if(input$useConsensus == "hit list") data <- outliers()
  else data <- consensusHitList()
  
  result <- isolate({
    result <- getTargets(data, group.miRNAs=input$group.miRNAs, group.miRNAs.threshold=input$group.miRNAs.threshold, databases=input$selectedTargetDBs)
    #if(input$excludeDBcol) result <- result[,setdiff(colnames(result), c("db_list"))]
    
    if(nrow(result) == 0) return (data.frame(error="No target has been found. reduce stringency or increase number of hits."))
    if(input$group.miRNAs)
    {
      data$Sample <- sub("mir", "miR", data$Sample)
      result$miRNA_list <- as.character(lapply(lapply(str_split(result$miRNA_list, "/"), function(x){
        temp.result <- ""
        for(miR in x){
          if(input$colorizeInTargetList){
            cat.color <- switch(unique(subset(data, Sample == miR)$category),
                                "included" = "yellow",
                                "suppressor" = "red",
                                "promotor" = "blue") 
            temp.result <- paste(temp.result, "<p style='color:", cat.color, "'>", miR, "</p>", sep="")
          }
          else{
            if(nchar(temp.result) == 0) temp.result <- miR
            else temp.result <- paste(temp.result, "/", miR, sep="")
          } 
        }
        return(temp.result)
      }), paste, collapse=""))
      
      #count promotors (blue) and suppressors (red), build difference
      result$promotor_vs_suppressor <- sapply(result$miRNA_list, function(x){
        blue <- gregexpr("blue", x)[[1]]
        red <- gregexpr("red", x)[[1]]
        if(blue[1] != -1 && red[1] != -1) return(length(blue) - length(red))
        else if(blue[1] == -1) return(-length(red))
        else return(length(blue))
      })
    }
    return(result)
  })                 
  return(result)
})

#indicator_matrix 
targets.indicator.matrix <- reactive({
  mirna.targets <- targetsForInteractionGraph()
  indicator.matrix <- as.data.frame.matrix(table(mirna.targets[,c("gene_id", "mature_miRNA")]))
  indicator.matrix[indicator.matrix > 1] <- 1
  return(indicator.matrix)
})