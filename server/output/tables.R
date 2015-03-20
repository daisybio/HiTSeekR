# Hit List #
formattedTable <- reactive({  
  data <- outliers()
  if(is.null(data)) return(NULL)
  
  if(!input$show.sem.in.hits)
    data <- data %>% dplyr::select(-ends_with("_sem")) 
  
  if("mature_from" %in% colnames(data))
    data <- data %>% dplyr::select(-one_of("mature_from", "mature_to", "evidence", "experiment", "similarity"))
  
  data[data$category %in% c("promotor"),"category"] <- "<div style='background:#80B1D3; text-align:center; border-radius: 15px; width:25px; height:25px;'>P</div>"
  data[data$category %in% c("suppressor"),"category"] <- "<div style='background:#FB8072; text-align:center; border-radius: 15px; width:25px; height:25px;'>S</div>"
  data[data$category %in% c("included"),"category"] <- "<div style='background:#FDB462; text-align:center; border-radius: 15px; width:25px; height:25px;'>I</div>"
  return(as.data.frame(data))
})

output$table_hits <- renderDataTable(formattedTable(), escape=FALSE)

# Raw data
output$table_rawData <- renderDataTable(rawData(), escape=FALSE)

# Processed data 
output$table_processedData <- renderDataTable(processedData(), escape=FALSE)

# Consensus hit list #
output$consensusHitList <- renderChart2({
  data <- consensusHitList()
  data[data$category %in% c("promotor"),"category"] <- "<div style='background:#80B1D3; text-align:center; border-radius: 15px; width:25px; height:25px;'>P</div>"
  data[data$category %in% c("suppressor"),"category"] <- "<div style='background:#FB8072; text-align:center; border-radius: 15px; width:25px; height:25px;'>S</div>"
  
  dTable(data, sPaginationType='full_numbers')
})

# miRNA targets #
output$mirna.targets.table <- renderDataTable(mirna.targets(), escape=FALSE)

# miRNA target permutation test results table
output$mirna.target.permutation.table <- renderDataTable(filtered.mirna.target.permutation(), escape=FALSE)

# Family hit rate # 
output$family.hitrate <- renderDataTable(family.hitrate(), escape=FALSE)

# GO enrichment analysis based on mRNA targets #
output$goEnrichmentTable <- renderChart2({
    dTable(goEnrichment(), sPaginationType='full_numbers')
})

# mirCancerDB #
output$mircancer.table <- renderDataTable({
  mirc.hits <- hits.mircancer()
  mirc.hits[mirc.hits$category %in% c("promotor"),"category"] <- "<div style='background:#80B1D3; text-align:center; border-radius: 15px; width:25px; height:25px;'>P</div>"
  mirc.hits[mirc.hits$category %in% c("suppressor"),"category"] <- "<div style='background:#FB8072; text-align:center; border-radius: 15px; width:25px; height:25px;'>S</div>"  
  mirc.hits$mirna_id <- paste("<a href='http://mircancer.ecu.edu/search.jsp?mirId=", mirc.hits$mirna_id, "&logic=&condition=And&cancerName=", mirc.hits$Cancer, "'>", mirc.hits$mirna_id, "</a>", sep="")
  return(mirc.hits)
}, escape=FALSE)
