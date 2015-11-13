output$uiOutput_drug_targets <- renderUI({
  if(is.null(input$screenType) && input$screenType != "compound"){
    stop("drug target prediction is only available for small compound / drug screens")
  }
  
  elements <- list(
    tabPanel("Drug target genes",  
             sidebarPanel(
               selectInput("selectedDrugTargetDBs", "currently selected:", c("STITCH"), "STITCH"),
               numericInput("drugTargetCutoff", "Cutoff for the STITCH combined score", value=500),
               helpText("STITCH Chemical-Protein Interactions database http://stitch.embl.de/")             

               #conditionalPanel(
               # condition = "input.selectedTargetDBs!='RNAhybrid_hsa'",
               # checkboxInput("get.gene.symbols", "Fetch gene symbols?", FALSE)
               #)
               #sliderInput("group.miRNAs.threshold", "List only genes that are targeted by x miRNAs", min=1, max=50, value=2, step=1),
               #sliderInput("at.least.hits", "List only targets found at least x times in databases", min = 1, max = 1000, value = 1, step = 1),
               #sliderInput("at.least", "List only targets found in at least x databases", min = 1, max = 6, value = 3, step = 1),
               #checkboxInput("excludeDBcol", "Exclude database text column from gene target list", TRUE),
               #actionButton("updateTargets", "Update miRNA target list")
             ),mainPanel( 
             dataTableOutput("drug.targets.table"),
             downloadButton('downloadDrugTargets', 'Download drug target list')
             )
    )#,
  #tabPanel("Interaction Graph", plotOutput("interactionGraph")),
  #tabPanel("miRcancer DB", shinyalert("mircancer_status"), dataTableOutput("mircancer.table"))
)

do.call(tabsetPanel, elements)
})


