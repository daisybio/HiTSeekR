output$uiOutput_consensus_hits_options <- renderUI({
  elements <- list(
      checkboxGroupInput("multiNormalizations", "Selection for consensus normalization:", normalizationChoices(), selected=c("rcentered", "rzscore", "Bscore")),
      sliderInput("multiThreshold", "Threshold for consensus normalization:", min = 1, max = 8, value= 3, step= 1)
  )
  do.call(sidebarPanel, elements)
})

output$uiOutput_consensus_hits <- renderUI({
  elements <- list(
    tabPanel("Consensus Hits List", 
             dataTableOutput("consensusHitList"), 
             downloadButton('downloadConsensusHits', 'Download')),    
    tabPanel("Consensus Hits Plot", showOutput("normcomparison", "polycharts")),
    tabPanel("Consensus Venn Diagram", plotOutput("consensusVennDiagram", height=400))  
  )
  do.call(tabsetPanel, elements)
})