output$interactionTable <- renderChart2({
  dTable(targetsForInteractionGraph())
})

# Hit List #
formattedTable <- reactive({  
  data <- outliers()
  
  if(is.null(data)) return(NULL)
  
  if(!input$show.sem.in.hits)
    data <- data %>% dplyr::select(-ends_with(".sem")) 
  
  data <- data %>% dplyr::select(-one_of("mature_from", "mature_to", "evidence", "experiment", "similarity"))
  data[data$category %in% c("promotor"),"category"] <- "<div style='background:#80B1D3; text-align:center; border-radius: 15px; width:25px; height:25px;'>P</div>"
  data[data$category %in% c("suppressor"),"category"] <- "<div style='background:#FB8072; text-align:center; border-radius: 15px; width:25px; height:25px;'>S</div>"
  data[data$category %in% c("included"),"category"] <- "<div style='background:#FDB462; text-align:center; border-radius: 15px; width:25px; height:25px;'>I</div>"
  return(data)
})

output$table_hits <- renderDataTable(formattedTable())

# Raw data
output$table_rawData <- renderDataTable(rawData())

# Processed data 
output$table_processedData <- renderDataTable(processedData())

# Consensus hit list #
output$consensusHitList <- renderChart2({
  data <- consensusHitList()
  data[data$category %in% c("promotor"),"category"] <- "<div style='background:#80B1D3; text-align:center; border-radius: 15px; width:25px; height:25px;'>P</div>"
  data[data$category %in% c("suppressor"),"category"] <- "<div style='background:#FB8072; text-align:center; border-radius: 15px; width:25px; height:25px;'>S</div>"
  
  dTable(data, sPaginationType='full_numbers')
})

# mRNA targets #
output$mirna.targets.table <- renderDataTable(targetsForInteractionGraph())

# Family hit rate # 
output$family.hitrate <- renderDataTable(family.hitrate())

# GO enrichment analysis based on mRNA targets #
output$goEnrichmentTable <- renderChart2({
  if(input$group.miRNAs && input$colorizeInTargetList){
    dTable(goEnrichment(), sPaginationType='full_numbers')
  } else{
    dTable(data.frame(error="This feature only works when miRNAs are grouped and colorized under 'Target Genes'!"))
  }
})

# mirCancerDB #
output$mircancerTable <- renderChart2({
  dTable(outliers.mircancer(), sPaginationType='full_numbers')
})
