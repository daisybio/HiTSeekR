output$uiOutput_mirna_targets <- renderUI({
  if(is.null(input$screenType) && input$screenType != "miRNA"){
    stop("miRNA target prediction is only available for miRNA inhibitor / mimics screens")
  }
  
  elements <- list(
    tabPanel("miRNA target genes",  
             wellPanel(
               selectInput("useConsensus", "Use normal hit list or consensus hit list for target identification?", c("hit list", "consensus hit list")),                    
               selectInput("selectedTargetDBs", "currently selected:", c(dbListTables(RmiR.Hs.miRNA_dbconn()), "RNAhybrid_hsa"), "RNAhybrid_hsa"),
               #helpText("tarbase is a database of experimentally verified targets. Other DBs deliver prediction based targets.")             
               conditionalPanel(
                 condition = "input.selectedTargetDBs=='RNAhybrid_hsa'",
                 sliderInput("rnah.p.value.threshold", "p-value threshold", min=-5, max=0, step=1, value=-3, ticks= c(0.00001, 0.0001, 0.001, 0.01, 0.1, 0))
               ),
               #sliderInput("group.miRNAs.threshold", "List only genes that are targeted by x miRNAs", min=1, max=50, value=2, step=1),
               #sliderInput("at.least.hits", "List only targets found at least x times in databases", min = 1, max = 1000, value = 1, step = 1),
               #sliderInput("at.least", "List only targets found in at least x databases", min = 1, max = 6, value = 3, step = 1),
               #checkboxInput("excludeDBcol", "Exclude database text column from gene target list", TRUE),
               actionButton("updateTargets", "Update miRNA target list")
             ), 
             dataTableOutput("mirna.targets.table"),
             downloadButton('downloadTargets', 'Download miRNA target list'),
             downloadButton('downloadHotnetGeneList', 'Download hotnet2 heat scores')
             
    ),
  tabPanel("miRNA target permutation test", 
    actionButton("mirna.target.permutation.button", "Start miRNA target permutation test"),
    numericInput("mirna.target.permutations", "Number of permutations", value=100, min=10, max=10000),    
    sliderInput("mirna.target.permutation.num.of.mirnas.cutoff", "Minimal number of miRNAs from hit list targeting a gene", value=1, min = 0, max=100, step=1),
    sliderInput("mirna.target.permutation.padj.cutoff", "adjusted p-value threshold", min=-5, max=0, step=1, value=-3, ticks= c(0.00001, 0.0001, 0.001, 0.01, 0.1, 0)),
    dataTableOutput("mirna.target.permutation.table"),
    downloadButton("downloadTargetPermutationTestResult", "Download")
  ),
  #tabPanel("Interaction Graph", plotOutput("interactionGraph")),
  tabPanel("Target Gene Enrichment", uiOutput("uiOutput_KPM"), mainPanel(
    shinyalert("kpm_status"),    
    checkboxInput("kpm_debug", "Show debug console", FALSE),
    conditionalPanel("input.kpm_debug", 
                     verbatimTextOutput("KPM.test")
    ),
    checkboxInput("kpm_d3", "Force directed network?", FALSE),
    conditionalPanel("input.kpm_d3",
                     sliderInput("highlight.kpm_d3", "Highlight genes in green that appear in more than x solution", min=1, max=20, value=5),
                     forceNetworkOutput("KPM.plot.d3")
    ),
    conditionalPanel("!input.kpm_d3",
                     plotOutput("KPM.plot.igraph", height=800, width=1200)
    )    
  )),
  tabPanel("miRcancer DB", shinyalert("mircancer_status"), dataTableOutput("mircancer.table"))
)

do.call(tabsetPanel, elements)
})

output$uiOutput_KPM <- renderUI({
  
  elements <- list(
    HTML('<img src="KPM_banner.png"/><br/><br/>'),
    textInput("kpm_URL", "KPM-Web URL:", "http://localhost:8080/kpm-web/"),  
    selectInput("kpm_strategy", "Strategy:", c("GLONE", "INES")),    
    selectInput("kpm_algorithm", "Algorithm:", list("Greedy"="Greedy", "Exact (FPT)"="Exact", "Ant Colony Optimization" = "ACO")),
    selectInput("kpm_network", "Network", KPM.network.list()),
    checkboxInput("kpm_ben_removal", "Remove border exception nodes?", TRUE),
    checkboxInput("random.miRNA.test", "Repeat with random miRNA set to estimate p-values", FALSE),
    numericInput("random.miRNA.iterations", "Iterations", min=10, max=1e6, value=10),
    checkboxInput("kpm_ranged", "ranged K and L?", FALSE),
    conditionalPanel("!input.kpm_ranged",
      conditionalPanel("input.kpm_strategy!='GLONE'",
        numericInput("kpm_K", "K (# node exceptions)", 1, min = 1, max = 100)        
      ),
      numericInput("kpm_L", "L (# case exceptions)", 1, min = 1, max = 1000)
    ),
    conditionalPanel("input.kpm_ranged",
      conditionalPanel("input.kpm_strategy!='GLONE'",
        numericInput("kpm_lower_K", "Lower limit for K", min=0, value=0),
        numericInput("kpm_upper_K", "Upper limit for K", min=0, value=5),
        numericInput("kpm_step_K", "Step size for K", min=1, value=1)
      ),      
      numericInput("kpm_lower_L", "Lower limit for L", min=0, value=0),
      numericInput("kpm_upper_L", "Upper limit for L", min=1, value=10),
      numericInput("kpm_step_L", "Step size for L", min=1, value=1)
    ),
    sliderInput("kpm_pathways", "Number of pathways", min = 1, max = 100, value=20),
    checkboxInput("kpm_perturbation", "Perturbe network?", FALSE),
    actionButton("startKPMButton", "Start KPM"), downloadButton('downloadIndicatorMatrix')
  )
  do.call(sidebarPanel, elements)
})


