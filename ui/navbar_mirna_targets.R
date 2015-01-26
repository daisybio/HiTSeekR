output$uiOutput_mirna_targets <- renderUI({
  #if(is.null(input$screenType) && input$screenType != "miRNA"){
  #  return(wellPanel("miRNA target prediction is only available for miRNA inhibitor / mimics screens"))
  #}
  
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
             downloadButton('downloadTargets', 'Download miRNA target list')
             
    ),
  #tabPanel("Interaction Graph", plotOutput("interactionGraph")),
  tabPanel("KeyPathwayMiner Target Network Enrichment", uiOutput("uiOutput_KPM"), mainPanel(
    shinyalert("kpm_status"),    
    checkboxInput("kpm_debug", "Show debug console", FALSE),
    conditionalPanel("input.kpm_debug", 
                     verbatimTextOutput("KPM.test")
    ),plotOutput("KPM.plot", height=800, width=1200))
  ),
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
    checkboxInput("kpm_ben_removal", "Remove border exception nodes?", TRUE),
    numericInput("kpm_K", "K (# node exceptions)", 1, min = 1, max = 100),
    numericInput("kpm_L", "L (# case exceptions)", 1, min = 1, max = 1000),
    actionButton("startKPMButton", "Start KPM"), downloadButton('downloadIndicatorMatrix')
  )
  do.call(sidebarPanel, elements)
})


