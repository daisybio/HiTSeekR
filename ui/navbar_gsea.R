output$uiOutput_htsanalyzer <- renderUI({
  if(input$startHTSanalyzer == 0) return(NULL)
  
  isolate({
    elements <- foreach(geneset.type = input$htsanalyzer.geneset.types) %do%
      tabPanel(geneset.type, 
               dataTableOutput(paste("htsanalyzer.results.table.", geneset.type, sep="")),
               downloadButton(paste("htsanalyzer.results.download.", geneset.type, sep=""))
      )
    do.call(tabsetPanel, elements)
  })
})

output$uiOutput_htsanalyzerOptions <- renderUI({
  
  elements <- list(
    if(input$screenType=="miRNA")
    {
      selectInput("htsanalyzer.miRNA.list", "Select miRNA target gene list", c("permutation test"="miRNA_permutation", "network enrichment"="miRNA_KPM"), "miRNA_permutation")
    } else
    {  
      selectInput("htsanalyzer.useConsensus", "Use consensus hit list for target identification?", c("hit list", "consensus hit list"))                    
    },
    selectInput("htsanalyzer.species", "Species:", 
           c("Drosophila melanogaster" = "Dm", 
             "Homo sapiens" = "Hs",
             "Rattus norvegicus" = "Rn",
             "Mus musculus" = "Mm",
             "Caenorhabditis elegans" = "Ce"
           ), selected = "Hs"         
    ),    
    checkboxGroupInput("htsanalyzer.geneset.types", "Select gene sets:", c("GO cellular compartment" = "GO_CC",
                                                                    "GO molecular function" = "GO_MF",
                                                                    "GO biological process" = "GO_BP",
                                                                    "KEGG pathways" ="PW_KEGG"),
                       selected = c("GO_CC", "GO_MF", "GO_BP", "PW_KEGG")
                       ),
    numericInput("htsanalyzer.pval.cutoff", "p-value cutoff", min = 0, max = 1, value= 0.05),
    numericInput("htsanalyzer.minimum.gene.set.size", "Minimal gene set size", min = 1, value = 10),    
    if(input$screenType=="miRNA")
    {
      HTML("Gene set enrichment analysis is not available for miRNA target lists")      
    } else
    {  
      checkboxInput("htsanalyzer.doGSEA", "Perform gene set enrichtment analysis", TRUE)
    },
    conditionalPanel(condition = "input['htsanalyzer.doGSEA']",
                     numericInput("htsanalyzer.permutations", "Number of permutations for GSEA", min=10, max=1000, value=100)
    ),          
    selectInput("htsanalyzer.adjust.method", "Method for p-value correction", p.adjust.methods, selected="BH"),
    actionButton("startHTSanalyzer", "Start Analysis")
  )
  do.call(sidebarPanel, elements)
})


