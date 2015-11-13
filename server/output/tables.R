# Hit List #
formattedTable <- function(exp.data, show.sd, show.all=TRUE, show.position=TRUE, rename=FALSE){  
  
  if(is.null(exp.data)) return(NULL)
  
  if(!show.sd){
    test.remove <- exp.data %>% dplyr::select(-ends_with("_sd")) 
    if(ncol(test.remove) > 0)
    {
      exp.data <- exp.data %>% dplyr::select(-ends_with("_sd")) 
    }
  }
  
  if(!show.all){
    exp.data <- exp.data %>% dplyr::select_(.dots = sapply(setdiff(normalizationChoices(), input$normalization), function(x){ paste("-", x, sep="")}))
  }
  
  if(!show.position){
    exp.data <- exp.data %>% dplyr::select(-Plate, -Well.position, -Row, -Column)
  }
  
  if("mature_from" %in% colnames(exp.data))
    exp.data <- exp.data %>% dplyr::select(-one_of("mature_from", "mature_to", "evidence", "experiment", "similarity"))
  
  if(nrow(exp.data) == 0) stop("No hits were found with the current settings")
  exp.data[exp.data$category %in% c("promotor"),"category"] <- "<div style='background:#80B1D3; text-align:center; border-radius: 15px; width:25px; height:25px;'>P</div>"
  exp.data[exp.data$category %in% c("suppressor"),"category"] <- "<div style='background:#FB8072; text-align:center; border-radius: 15px; width:25px; height:25px;'>S</div>"
  exp.data[exp.data$category %in% c("included"),"category"] <- "<div style='background:#FDB462; text-align:center; border-radius: 15px; width:25px; height:25px;'>I</div>"
    
  if(rename && input$screenType == "siRNA"){
    exp.data <- exp.data %>% dplyr::rename(Category = category, `Entrez ID` = gene_id, `Gene Symbol` = gene_symbol)
  }
  else
  {
    exp.data <- exp.data %>% dplyr::rename(Category = category)
  }
  return(as.data.frame(exp.data))
}

output$table_hits <- renderDataTable(formattedTable(outliers(), input$show.sd.in.hits, input$showAllScores, input$showSamplePosition, rename=TRUE), escape=FALSE)

# Legend for category
output$cat_legend <- renderText("<b>Categories:</b><br/><span>Promotor: <div style='background:#80B1D3; text-align:center; border-radius: 15px; width:25px; height:25px;'>P</div></span> 
<span>Suppressor: <div style='background:#FB8072; text-align:center; border-radius: 15px; width:25px; height:25px;'>S</div></span>
<span>Included by regular expression: <div style='background:#FDB462; text-align:center; border-radius: 15px; width:25px; height:25px;'>I</div></span>"
)
# Raw data
output$table_rawData <- renderDataTable(rawData(), escape=FALSE)

# Processed data 
output$table_processedData <- renderDataTable(processedData(), escape=FALSE)

# Consensus hit list #
#output$consensusHitList <- renderDataTable(formattedTable(consensusHitList(), input$show.sd.in.hits), escape=FALSE)

# miRNA targets #
output$mirna.targets.table <- renderDataTable(mirna.targets(), escape=FALSE)

# Drug targets #
output$drug.targets.table <- renderDataTable(drug.targets(), escape=FALSE)

# miRNA target permutation test results table
output$mirna.target.permutation.table <- renderDataTable({
  if(is.null(filtered.mirna.target.permutation())) return(NULL)
  filtered.mirna.target.permutation() %>% dplyr::select(-mature_miRNA)
  }, escape=FALSE)

# Family hit rate # 
output$family.hitrate <- renderDataTable(family.hitrate(), escape=FALSE)

# HTSAnalyzer results #
output$htsanalyzer.results.table.GO_CC <- renderDataTable(htsanalyzer.results()[["GO_CC"]], escape=FALSE)
output$htsanalyzer.results.table.GO_MF <- renderDataTable(htsanalyzer.results()[["GO_MF"]], escape=FALSE)
output$htsanalyzer.results.table.GO_BP <- renderDataTable(htsanalyzer.results()[["GO_BP"]], escape=FALSE)
output$htsanalyzer.results.table.PW_KEGG <- renderDataTable(htsanalyzer.results()[["PW_KEGG"]], escape=FALSE)
output$htsanalyzer.results.table.REACTOME <- renderDataTable(htsanalyzer.results()[["REACTOME"]], escape=FALSE)

# KPM selected graph - nodes
output$show_kpm_nodes <- renderDataTable( formattedTable(kpm.node.table(), TRUE), escape=FALSE)

# KPM selected gene info
output$kpm_gene_details <- reactive({
  table.data <- kpm.node.table()
  sel.gene <- table.data[input$show_kpm_nodes_row_last_clicked, "gene_id"]
  return(mygene::query(sel.gene, fields=c("go")))
})

output$kpm_gene_details_cc <- renderDataTable({
  kpm_gene_details()$hits$go$CC[[1]]
})

output$kpm_gene_details_bp <- renderDataTable({
  kpm_gene_details()$hits$go$BP[[1]]
})

output$kpm_gene_details_mf <- renderDataTable({
  kpm_gene_details()$hits$go$MF[[1]]
})

# mirCancerDB #
output$mircancer.table <- renderDataTable({
  mirc.hits <- hits.mircancer()
  mirc.hits[mirc.hits$category %in% c("promotor"),"category"] <- "<div style='background:#80B1D3; text-align:center; border-radius: 15px; width:25px; height:25px;'>P</div>"
  mirc.hits[mirc.hits$category %in% c("suppressor"),"category"] <- "<div style='background:#FB8072; text-align:center; border-radius: 15px; width:25px; height:25px;'>S</div>"  
  mirc.hits$mirna_id <- paste("<a target='_blank' href='http://mircancer.ecu.edu/search.jsp?mirId=", mirc.hits$mirna_id, "&logic=&condition=And&cancerName=", mirc.hits$Cancer, "'>", mirc.hits$mirna_id, "</a>", sep="")
  return(mirc.hits)
}, escape=FALSE)

# DIANA mirPATH table #
output$mirpath.table <- renderDataTable( mirpath.results() )
