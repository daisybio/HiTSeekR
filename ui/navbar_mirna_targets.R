output$uiOutput_mirna_targets <- renderUI({
  
  if(input$screenType != "miRNA"){
    return(wellPanel("miRNA target prediction is only available for miRNA inhibitor / mimics screens"))
  }

  elements <- list(
    tabPanel("Options", 
             selectInput("useConsensus", "Use normal hit list or consensus hit list for target identification?", c("hit list", "consensus hit list")),
             checkboxInput("showTargetDBs", "Select miRNA target databases", FALSE),
             conditionalPanel(
               condition = "input.showTargetDBs",
               checkboxGroupInput("selectedTargetDBs", "currently selected:", c(dbListTables(RmiR.Hs.miRNA_dbconn()), "RNAhybrid_hsa"), "RNAhybrid_hsa"),
               helpText("tarbase is a database of experimentally verified targets. Other DBs deliver prediction based targets.")
             ),
             checkboxInput("group.miRNAs", "Group miRNAs to gene targets", FALSE),
             checkboxInput("remove.nas.from.target.list", "Remove targets with NA entry as gene symbol", TRUE),                      
             checkboxInput("colorizeInTargetList", "Colorize and count suppressors vs promotors", FALSE),
             sliderInput("group.miRNAs.threshold", "List only genes that are targeted by x miRNAs", min=1, max=50, value=2, step=1),
             sliderInput("at.least.hits", "List only targets found at least x times in databases", min = 1, max = 1000, value = 1, step = 1),
             sliderInput("at.least", "List only targets found in at least x databases", min = 1, max = 6, value = 3, step = 1),
             checkboxInput("excludeDBcol", "Exclude database text column from gene target list", TRUE),
             actionButton("updateTargets", "Update miRNA target list")
    ),
    tabPanel("miRNA target genes", 
             dataTableOutput("mirna.targets.table"),
             downloadButton('downloadTargets', 'Download miRNA target list')
    ),
    tabPanel("Interaction Graph", plotOutput("interactionGraph"))            
  )

do.call(tabsetPanel, elements)
})