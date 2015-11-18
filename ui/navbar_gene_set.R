gene.sets <- c("GO cellular compartment" = "GO_CC",
               "GO molecular function" = "GO_MF",
               "GO biological process" = "GO_BP",
               "KEGG pathways" ="PW_KEGG",
               "REACTOME pathways" = "REACTOME")

output$uiOutput_htsanalyzer <- renderUI({
  if(input$startHTSanalyzer == 0) return(NULL)
  
  isolate({
    elements <- foreach(geneset.type = input$htsanalyzer.geneset.types) %do%
      tabPanel(names(gene.sets)[gene.sets == geneset.type], 
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
       selectInput("htsanalyzer.miRNA.list", "Select miRNA target gene list", c("miRNA target gene list"="miRNA_targets", "high confidence target genes"="miRNA_permutation"), "miRNA_permutation")#, "network enrichment genes"="miRNA_KPM"), "miRNA_permutation")
     }, 
#     selectInput("htsanalyzer.species", "Species:", 
#                 c("Drosophila melanogaster" = "Dm", 
#                   "Homo sapiens" = "Hs",
#                   "Rattus norvegicus" = "Rn",
#                   "Mus musculus" = "Mm",
#                   "Caenorhabditis elegans" = "Ce"
#                 ), selected = "Hs"         
#     ),    
    checkboxGroupInput("htsanalyzer.geneset.types", "Select gene sets:", gene.sets,
                       selected = c("GO_CC", "GO_MF", "GO_BP", "PW_KEGG", "REACTOME")
    ),
    numericInput("htsanalyzer.pval.cutoff", "Adjusted p-value cutoff", min = 0, max = 1, value= 0.05),
    numericInput("htsanalyzer.minimum.gene.set.size", "Minimal gene set size", min = 1, value = 10),    
    if(input$screenType %in% c("miRNA", "compound"))
    {
      HTML("Note: gene set enrichment analysis is not available for miRNA or compound target lists.")      
    } else 
    {  
      checkboxInput("htsanalyzer.doGSEA", "Perform gene set enrichtment analysis", FALSE)
    },
    conditionalPanel(condition = "input['htsanalyzer.doGSEA']",
                     numericInput("htsanalyzer.permutations", "Number of permutations for GSEA", min=10, max=gsea.max.permutations, value=100)
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
               HTML('<div class="shinyalert alert fade alert-info in">Gene set analysis is perfomed based on pre-defind gene sets. The enrichment of hit genes in a particular gene set can be tested for significance using hyper-geometric tests or gene set enrichment analysis.</div>'), 
               shinyalert("htsanalyzer_status"),
               conditionalPanel("input['htsanalyzer.doGSEA']",
               selectInput("htsanalyzer.resultType", "Select results", 
                           c("Hypergeometric Test"= "HyperGeo.results", 
                             "Gene set enrichment analysis" = "GSEA.results",
                             "Significant adjusted p-values in both" = "Sig.adj.pvals.in.both")
               )),
               uiOutput("uiOutput_htsanalyzer")
             )
    ),
    tabPanel("De novo network enrichment", uiOutput("uiOutput_KPM"), mainPanel(
      HTML('<div class="shinyalert alert fade alert-info in">In contrast to gene set analysis, de novo network enrichment analysis does not rely on pre-defined gene sets, but directly extracts enriched sub-networks from large gene / protein interaction networks.</div>'), 
      shinyalert("kpm_status"),    
      #textOutput("ind.matrix.props"),    
      #checkboxInput("kpm_debug", "Show debug console", FALSE),
      #conditionalPanel("input.kpm_debug", 
      #                 verbatimTextOutput("KPM.test")
      #),
      #checkboxInput("kpm_d3", "Force directed network?", FALSE),
      #conditionalPanel("input.kpm_d3",
      #                 sliderInput("highlight.kpm_d3", "Highlight genes in green that appear in more than x solution", min=1, max=20, value=5),
      #                 forceNetworkOutput("KPM.plot.d3")
      #),
      #conditionalPanel("!input.kpm_d3",
      uiOutput("uiOutput_kpm_plot")
      #)    
    ))
  )  
  if(input$screenType == "miRNA") 
    elements <- c(elements, list(tabPanel("DIANA mirPATH", 
                                          sidebarPanel(
                                            HTML('<div class="shinyalert alert fade alert-info in">For details regarding DIANA and the settings shown below click <a href=\'http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=mirpath\' target=\'_blank\'><u>here</u></a>.</div>'), 
                                            selectInput("mirpath_selection", "Gene selection method", c("Genes union" = 0, "Genes intersection" = 1, "Pathways union" = 2, "Pathways intersection" = 3), 2),
                                            conditionalPanel("input.mirpath_selection == 1", 
                                                             numericInput("mirpath_cutoff", "Cutoff for gene intersection", min = 0, value = 0, max=100)
                                            ),
                                            numericInput("mirpath_fdr", "Multiple testing (FDR) cutoff", min=0, max=1, step=0.01, value=0.05),
                                            checkboxInput("mirpath_conservative", "Use conservative statistics", value=TRUE),
                                            #selectInput("mirpath_method", "miRNA target prediction method", c("tarbase", "microT-CDS"), value="microT-CDS"),
                                            numericInput("mirpath_threshold", "microT-CDS score cutoff", min=0.7, max=1.0, value=0.8, step=0.01)
                                          ),
                                          mainPanel(dataTableOutput("mirpath.table"), downloadButton("downloadMirPathResults"))
                                          )
                                 )
                  )
  
  if(input$isHitList){
    elements[[1]] <-  tabPanel("Gene set analysis",                               
                                    HTML('<div class="shinyalert alert fade alert-info in">Gene set analysis cannot be performed directly on a hit list, since information about the entire experiment is needed to calculate the necessary p-values.</div>'))
                               
  }
  do.call(tabsetPanel, elements)
})


#KPM main panel with result plots
output$uiOutput_kpm_plot <- renderUI({
  kpm.res <- KPM.result()
  if(is.null(kpm.res)) return(NULL)
  
  elements <- list(
    wellPanel(
      checkboxInput("kpm_union_graph", "Show union graph of top 20 extracted key pathways?", value = FALSE),
      conditionalPanel(condition = '!input.kpm_union_graph',
        sliderInput("kpm_selected_solution", "Select which of the top 20 key pathways to display", min=1, max=20, value=1)
      )
    ),
    tabsetPanel(
      tabPanel("Graph",
        plotOutput("KPM.plot.igraph", height=800, width=900),
        downloadButton('download_kpm_SIF', 'Download this graph as Cytoscape compatible SIF file')
      ),  
      tabPanel("Table", 
       dataTableOutput("show_kpm_nodes"),
       downloadButton('download_kpm_node_table', 'Download')
       #tabsetPanel(position = "left", tabPanel("GO Biological Process",
       #dataTableOutput("kpm_gene_details")))
      )
    )  
  )
  do.call(mainPanel, elements)
})

#KPM sidebar panel with options
output$uiOutput_KPM <- renderUI({
  elements <- list(
    HTML('<img src="KPM_banner.png"/><br/><br/>'),
    #textInput("kpm_URL", "KPM-Web URL:", "http://localhost:8080/kpm-web/"),  
    HTML('<div class="shinyalert alert fade alert-info in">For details regarding KeyPathwayMiner and the settings shown below click <a href=\'http://tomcat.compbio.sdu.dk/keypathwayminer/\' target=\'_blank\'><u>here</u></a>.</div>'), 
    if(input$screenType=="miRNA")
    {
      selectInput("KPM.miRNA.list", "Select miRNA target gene list", c("miRNA target gene list"="miRNA_targets", "high confidence target genes"="miRNA_permutation"), "miRNA_targets")
    },
    #selectInput("kpm_strategy", "Strategy:", c("GLONE", "INES"), value="INES"),    
    #selectInput("kpm_algorithm", "Algorithm:", list("Greedy"="Greedy", "Exact (FPT)"="Exact", "Ant Colony Optimization" = "ACO")),
    selectInput("kpm_network", "Network", KPM.network.list()),
    #conditionalPanel(condition="input.kpm_strategy=='INES'",
      checkboxInput("kpm_ben_removal", "Remove border exception nodes?", TRUE),
    #),    
    #checkboxInput("kpm_ranged", "ranged K and L?", FALSE),
    #conditionalPanel("!input.kpm_ranged",
    if(input$screenType == "miRNA"){
      sliderInput("kpm_L", paste("How many times should a gene be targeted at least to be considered?"), 1, min = 1, max = length(unique(mirna.targets()$mature_miRNA)), step=1)
    }
    else if(input$screenType == "compound")
    {
      sliderInput("kpm_L", paste("How many times should a gene be targeted at least to be considered?"), 1, min = 1, max = length(unique(drug.targets()$PubChem_CID)), step=1)
    },
    #conditionalPanel(condition="input.kpm_strategy!='GLONE'",
    sliderInput("kpm_K", paste("Number of exception genes allowed"), 1, min = 0, max = 3, step=1),        
    #),
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
    #sliderInput("kpm_pathways", "Number of pathways", min = 1, max = 100, value=20),
    #checkboxInput("kpm_perturbation", "Perturbe network?", FALSE),
    actionButton("startKPMButton", "Start KeyPathwayMiner Analysis", styleclass="primary"), 
    hr(),
    HTML('<div class="shinyalert alert fade alert-info in">Alternatively, you can also download the indicator matrix needed as input for KeyPathwayMiner and perform the analysis directly in <a href=\'http://tomcat.compbio.sdu.dk/keypathwayminer/\' target=\'_blank\'><u>KeyPathwayMiner Web</u></a> or using the <a href="http://apps.cytoscape.org/apps/keypathwayminer" target="_blank"><u>Cytoscape app</u></a>.</div>'), 
    downloadButton('downloadIndicatorMatrix', 'Download indicator matrix')
  )
  do.call(sidebarPanel, elements)
})
