output$uiOutput_geneOntology <- renderUI({

elements <- list(
  sliderInput("goSignThreshold", "Select promotor_vs_suppressor threshold for significant genes:", min = 1, max = 15, value = 3, step = 1),
  sliderInput("goTopNodes", "List top x GO terms", min=10, value=100, max=1000, step=10),
  selectInput("goUpOrDown", "Promotors or Suppressors", choices = c("promotors", "suppressors", "both")),
  selectInput("goDomain", "Select ontology:", choices = list("biological process" = "BP", "molecular function" = "MF", "cellular location" = "CL") ),
  selectInput("goOrderMethod", "Order result by which method?", choices = list("Kolmogorov-Smirnov (elim)" = "elim.KS", "Kolmogorov-Smirnov (classic)" = "KS", "Fisher" = "Fisher")),
  selectInput("goSelectedMethod", "Select method for the graph:", choices=c("elim")),
  selectInput("goUseInfo", "Select info to be displayed in the graph:", choices= c("all", "def")),
  sliderInput("goSelectedNodes", "Select number of GO terms in the graph", min=1, max=25, step=1, value=5),
  chartOutput("goEnrichmentTable", "datatables"),
  downloadButton('dlGnOntGraph', 'Download GO enrichment graph as PDF'),
  downloadButton('dlGnOntTbl', 'Download GO enrichment table')
)
do.call(tabPanel, elements)
})