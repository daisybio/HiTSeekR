output$uiOutput_consensus_hits <- renderUI({
  elements <- list(
    tabPanel("Consensus Hits Plot",
           sidebarPanel(
             checkboxGroupInput("multiNormalizations", "Selection for consensus normalization:", c("Raw (CTB)" = "CTB", "Percentage of Control" = "poc", "Normalized percentage inhibition" = "npi", "Centered with mean" = "centered", "Centered with median" = "rcentered", "z-score" = "zscore", "robust z-score" = "rzscore", "B-score"="Bscore", "Position effect normalization"="posEffectNorm"), selected=c("rcentered", "rzscore", "Bscore", "posEffectNorm")),
             sliderInput("multiThreshold", "Threshold for consensus normalization:", min = 1, max = 8, value= 3, step= 1)
           ),
           showOutput("normcomparison", "polycharts"),
           downloadButton('downloadConsensusHits', 'Download consensus hit list')
    ),
    tabPanel("Consensus Hits List", chartOutput("consensusHitList", "datatables")),
    tabPanel("Consensus Venn Diagram", plotOutput("consensusVennDiagram", height=400))  
  )
  do.call(tabsetPanel, elements)
})