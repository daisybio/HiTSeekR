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
      selectInput("htsanalyzer.miRNA.list", "Select miRNA target gene list", c("miRNA target gene list"="miRNA_targets", "high confidence target genes"="miRNA_permutation", "network enrichment genes"="miRNA_KPM"), "miRNA_permutation")
    } else
    {  
      selectInput("htsanalyzer.useConsensus", "Use consensus hit list for target identification?", c("hit list", "consensus hit list"))                    
    },
#     selectInput("htsanalyzer.species", "Species:", 
#                 c("Drosophila melanogaster" = "Dm", 
#                   "Homo sapiens" = "Hs",
#                   "Rattus norvegicus" = "Rn",
#                   "Mus musculus" = "Mm",
#                   "Caenorhabditis elegans" = "Ce"
#                 ), selected = "Hs"         
#     ),    
    checkboxGroupInput("htsanalyzer.geneset.types", "Select gene sets:", c("GO cellular compartment" = "GO_CC",
                                                                           "GO molecular function" = "GO_MF",
                                                                           "GO biological process" = "GO_BP",
                                                                           "KEGG pathways" ="PW_KEGG"),
                       selected = c("GO_CC", "GO_MF", "GO_BP", "PW_KEGG")
    ),
    numericInput("htsanalyzer.pval.cutoff", "p-value cutoff", min = 0, max = 1, value= 0.05),
    numericInput("htsanalyzer.minimum.gene.set.size", "Minimal gene set size", min = 1, value = 10),    
    if(input$screenType %in% c("miRNA", "compound"))
    {
      HTML("Note: gene set enrichment analysis is not available for miRNA or compound target lists.")      
    } else 
    {  
      checkboxInput("htsanalyzer.doGSEA", "Perform gene set enrichtment analysis", TRUE)
    },
    conditionalPanel(condition = "input['htsanalyzer.doGSEA']",
                     numericInput("htsanalyzer.permutations", "Number of permutations for GSEA", min=10, max=1000, value=100)
    ),          
    #selectInput("htsanalyzer.adjust.method", "Method for p-value correction", p.adjust.methods, selected="BH"),
    actionButton("startHTSanalyzer", "Start Analysis", styleclass="primary")
  )
  do.call(sidebarPanel, elements)
})



output$uiOutput_gene_set_analysis <- renderUI({    
  elements <- list(
    tabPanel("Gene set analysis", uiOutput("uiOutput_htsanalyzerOptions"), 
             mainPanel(
               selectInput("htsanalyzer.resultType", "Select results", 
                           c("Hypergeometric Test"= "HyperGeo.results", 
                             "Gene set enrichment analysis" = "GSEA.results",
                             "Significant p-values in both" = "Sig.pvals.in.both",
                             "Significant adjusted p-values in both" = "Sig.adj.pvals.in.both")
               ), 
               uiOutput("uiOutput_htsanalyzer")
             )
    ),
    tabPanel("De novo network enrichment", uiOutput("uiOutput_KPM"), mainPanel(    
      shinyalert("kpm_status"),    
      #textOutput("ind.matrix.props"),    
      #checkboxInput("kpm_debug", "Show debug console", FALSE),
      #conditionalPanel("input.kpm_debug", 
      #                 verbatimTextOutput("KPM.test")
      #),
      checkboxInput("kpm_d3", "Force directed network?", FALSE),
      conditionalPanel("input.kpm_d3",
                       sliderInput("highlight.kpm_d3", "Highlight genes in green that appear in more than x solution", min=1, max=20, value=5),
                       forceNetworkOutput("KPM.plot.d3")
      ),
      conditionalPanel("!input.kpm_d3",
                       plotOutput("KPM.plot.igraph", height=800, width=1200)
      )    
    ))
  )  
    do.call(tabsetPanel, elements)
})


output$uiOutput_KPM <- renderUI({
  
  elements <- list(
    HTML('<img src="KPM_banner.png"/><br/><br/>'),
    #textInput("kpm_URL", "KPM-Web URL:", "http://localhost:8080/kpm-web/"),  
    if(input$screenType=="miRNA")
    {
      selectInput("KPM.miRNA.list", "Select miRNA target gene list", c("miRNA target gene list"="miRNA_targets", "high confidence target genes"="miRNA_permutation", "network enrichment genes"="miRNA_KPM"), "miRNA_permutation")
    } else
    {  
      selectInput("KPM.useConsensus", "Use consensus hit list for target identification?", c("hit list", "consensus hit list"))                    
    },
    selectInput("kpm_strategy", "Strategy:", c("GLONE", "INES")),    
    selectInput("kpm_algorithm", "Algorithm:", list("Greedy"="Greedy")), #"Exact (FPT)"="Exact", "Ant Colony Optimization" = "ACO")),
    selectInput("kpm_network", "Network", KPM.network.list()),
    conditionalPanel(condition="input.kpm_strategy=='INES'",
      checkboxInput("kpm_ben_removal", "Remove border exception nodes?", TRUE)
    ),    
    #checkboxInput("kpm_ranged", "ranged K and L?", FALSE),
    #conditionalPanel("!input.kpm_ranged",
    conditionalPanel(condition="input.kpm_strategy!='GLONE'",
         numericInput("kpm_K", paste("K (# node exceptions)"), 1, min = 1, max = 100)        
    ),
    numericInput("kpm_L", paste("L (# case exceptions)"), 1, min = 1, max = 1000),
    #),
    #conditionalPanel("input.kpm_ranged",
    #  conditionalPanel("input.kpm_strategy!='GLONE'",
    #    numericInput("kpm_lower_K", "Lower limit for K", min=0, value=0),
    #    numericInput("kpm_upper_K", "Upper limit for K", min=0, value=5),
    #    numericInput("kpm_step_K", "Step size for K", min=1, value=1)
    #  ),      
    #  numericInput("kpm_lower_L", "Lower limit for L", min=0, value=0),
    #  numericInput("kpm_upper_L", "Upper limit for L", min=1, value=10),
    #  numericInput("kpm_step_L", "Step size for L", min=1, value=1)
    #),
    sliderInput("kpm_pathways", "Number of pathways", min = 1, max = 100, value=20),
    #checkboxInput("kpm_perturbation", "Perturbe network?", FALSE),
    actionButton("startKPMButton", "Start KPM"), downloadButton('downloadIndicatorMatrix')
  )
  do.call(sidebarPanel, elements)
})
