## find mRNA targets ##
mirna.targets <- reactive({
  
  input$updateTargets
  
  data <- selectedHitList()
  
  result <- isolate({
    
    #check if accession type is mirna_id we need to convert to mature mirna ids in that case
    if(input$accessionColType == "MI")
    {
      showshinyalert(session, "mirna_target_status", "For pri-miRNA (Accession MI) it is unclear which of the two strands (3p or 5p) is active. Therefore, targets for both strands are predicted.","warning")
      mimat <- as.data.frame(mirbaseMATURE)
      data <- left_join(data, mimat, by=c("mirna_id"))
    }
    
    if(is.null(input$rnah.p.value.threshold)) rnah.p.value.threshold <- NULL
    else rnah.p.value.threshold <- 10^input$rnah.p.value.threshold
    
    result <- getTargets(data, rnah.pvalue.threshold=rnah.p.value.threshold, databases=input$selectedTargetDBs)
    #if(input$excludeDBcol) result <- result[,setdiff(colnames(result), c("db_list"))]
    
    if(nrow(result) == 0) stop("No target has been found. reduce stringency or increase number of hits.")
    
    return(result)
  })                 
  return(result)
})

#indicator_matrix 
targets.indicator.matrix <- reactive({
  mirna.targets <- mirna.targets()
  indicator.matrix <- as.data.frame.matrix(table(mirna.targets[,c("gene_id", "mature_miRNA")]))
  indicator.matrix[indicator.matrix > 1] <- 1
  return(indicator.matrix)
})

# mirna.targets.formatted <- reactive({
#   mirna.targets <- mirna.targets()
#   
#   if(group.miRNAs){
#     miRNA.targets <- miRNA.targets %>% 
#       group_by(mature_miRNA, gene_id) %>% 
#       summarize(count=n(), source_db=db) %>% 
#       filter(count > group.miRNAs.threshold)
#   }
#   #formatted table with promotors and suppressors and grouped genes
#   if(input$group.miRNAs)
#   {
#     data$Sample <- sub("mir", "miR", data$Sample)
#     result$miRNA_list <- as.character(lapply(lapply(str_split(result$miRNA_list, "/"), function(x){
#       temp.result <- ""
#       for(miR in x){
#         if(input$colorizeInTargetList){
#           cat.color <- switch(unique(subset(data, Sample == miR)$category),
#                               "included" = "yellow",
#                               "suppressor" = "red",
#                               "promotor" = "blue") 
#           temp.result <- paste(temp.result, "<p style='color:", cat.color, "'>", miR, "</p>", sep="")
#         }
#         else{
#           if(nchar(temp.result) == 0) temp.result <- miR
#           else temp.result <- paste(temp.result, "/", miR, sep="")
#         } 
#       }
#       return(temp.result)
#     }), paste, collapse=""))
#     
#     #count promotors (blue) and suppressors (red), build difference
#     result$promotor_vs_suppressor <- sapply(result$miRNA_list, function(x){
#       blue <- gregexpr("blue", x)[[1]]
#       red <- gregexpr("red", x)[[1]]
#       if(blue[1] != -1 && red[1] != -1) return(length(blue) - length(red))
#       else if(blue[1] == -1) return(-length(red))
#       else return(length(blue))
#     })
# })
#   