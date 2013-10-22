library(shiny)
require(rCharts)
require(ggplot2)
require(grid)
require(gridExtra)
require(scales)
require(RmiR)

shinyUI(pageWithSidebar(
  headerPanel("miRNAice"),
  sidebarPanel(
    selectInput("dataset", "Select a dataset", choices = c("MCF7wt miRNA inhibitor screen","MCF12Awt miRNA inhibitor screen", "own dataset")),
    conditionalPanel(
      condition="input.dataset=='own dataset'",
      fileInput("file", "Upload a pre-processed data set", multiple=FALSE)
    ),
    selectInput("normalization", "Raw Data / Normalization:", 
                choices = c("CTB", "poc", "npi", "centered", "rcentered", "zscore", "rzscore", "Bscore", "posEffectNorm")),
    conditionalPanel(
      condition = "input.normalization=='poc' || input.normalization=='npi'",
      helpText("Warning: When using percentage of control or normalized percentage inhibition: plate 12 is normalized using controls from plate 11!")
    ),
    selectInput("method", "Method for margin calculation:", choices = c("SD", "MAD", "quartile")),
    sliderInput("margin", 
                "Margin:", 
                min = 0.5,
                max = 5, 
                value = 2.5, step= 0.5),
    actionButton("updateNormalization", "Update Settings"),
    checkboxInput("includeControls", "Include control column in hit list", FALSE),
    checkboxInput("showSEM", "Show SEM in bubble charts", TRUE),
    checkboxInput("showFilterOptions", "Show filter options", FALSE),
    conditionalPanel(
      condition = "input.showFilterOptions",
      helpText("Filter options:"),
      textInput("exclude", "Exclude miRNAs by regular expression:", value=""),
      actionButton("updateExclusion", "Update exclusion filter"),
      textInput("include", "Include miRNAs in hit list by regular expression:", value=""),
      actionButton("updateInclusion", "Update inclusion filter"),
      helpText("For example, you can select all let-7 like this: let-7, or you can select several miRs like this: mir-(765|558)")
    ),
    helpText("Additional options:"),
    conditionalPanel(
      condition = "input.tabs1=='Heatmap'",
      selectInput("colorA", "Select low signal color for heatmap:", 
                choices = c("red", "blue", "darkblue", "steelblue", "magenta", "yellow", "white", "green")),
      selectInput("colorB", "Select high signal color for heatmap:", c("yellow", "red", "darkblue", "blue", "steelblue", "magenta", "white", "green"))
    ),
    conditionalPanel(
      condition = "input.tabs1=='Plate Viewer'",
      sliderInput("plateSelected", "Select a plate for the plate view:", min=1, max=12, value=1, step=1)
    ),
    conditionalPanel(
      condition = "input.tabs1=='Consensus Hits Plot' || input.tabs1=='Consensus Hits List' || input.tabs1=='Consensus Venn Diagram'",
      checkboxGroupInput("multiNormalizations", "Selection for consensus normalization:", c("Raw (CTB)" = "CTB", "Percentage of Control" = "poc", "Normalized percentage inhibition" = "npi", "Centered with mean" = "centered", "Centered with median" = "rcentered", "z-score" = "zscore", "robust z-score" = "rzscore", "B-score"="Bscore", "Position effect normalization"="posEffectNorm"), selected=c("Centered with median", "robust z-score", "B-score", "Position effect normalization")),
      sliderInput("multiThreshold", "Threshold for consensus normalization:", min = 1, max = 8, value= 3, step= 1)
    ),
    conditionalPanel(
      condition = "input.tabs1=='GO Enrichment'",
      sliderInput("goSignThreshold", "Select promotor_vs_suppressor threshold for significant genes:", min = 1, max = 15, value = 3, step = 1),
      sliderInput("goTopNodes", "List top x GO terms", min=10, value=100, max=1000, step=10),
      selectInput("goUpOrDown", "Promotors or Suppressors", choices = c("promotors", "suppressors", "both")),
      selectInput("goDomain", "Select ontology:", choices = list("biological process" = "BP", "molecular function" = "MF", "cellular location" = "CL") ),
      selectInput("goOrderMethod", "Order result by which method?", choices = list("Kolmogorov-Smirnov (elim)" = "elim.KS", "Kolmogorov-Smirnov (classic)" = "KS", "Fisher" = "Fisher")),
      selectInput("goSelectedMethod", "Select method for the graph:", choices=c("elim")),
      selectInput("goUseInfo", "Select info to be displayed in the graph:", choices= c("all", "def")),
      sliderInput("goSelectedNodes", "Select number of GO terms in the graph", min=1, max=25, step=1, value=5),
      downloadButton('dlGnOntGraph', 'Download GO enrichment graph as PDF'),
      downloadButton('dlGnOntTbl', 'Download GO enrichment table')
    ),
    conditionalPanel(
      condition = "input.tabs1=='Hits List'",
      downloadButton('downloadHits', 'Download hit list')
    ),  
    conditionalPanel(
      condition = "input.tabs1=='Consensus Hits List'",
      downloadButton('downloadConsensusHits', 'Download consensus hit list')
    ), 
    conditionalPanel(
      condition = "input.tabs1=='Target Genes'||input.tabs1=='Interaction Graph'||input.tabs1=='Interaction Table'",
      selectInput("useConsensus", "Use normal hit list or consensus hit list for target identification?", c("hit list", "consensus hit list")),
      checkboxInput("showTargetDBs", "Select miRNA target databases", FALSE),
      conditionalPanel(
        condition = "input.showTargetDBs",
        checkboxGroupInput("selectedTargetDBs", "currently selected:", dbListTables(RmiR.Hs.miRNA_dbconn()), dbListTables(RmiR.Hs.miRNA_dbconn())),
        helpText("tarbase is a database of experimentally verified targets. Other DBs deliver prediction based targets.")
      ),
      checkboxInput("group.miRNAs", "Group miRNAs to gene targets", FALSE),
      checkboxInput("colorizeInTargetList", "Colorize and count suppressors vs promotors", FALSE),
      sliderInput("group.miRNAs.threshold", "List only targets that are targeted by x miRNAs", min=1, max=50, value=2, step=1),
      sliderInput("at.least.hits", "List only targets found at least x times in databases", min = 1, max = 1000, value = 1, step = 1),
      sliderInput("at.least", "List only targets found in at least x databases", min = 1, max = 6, value = 3, step = 1),
      checkboxInput("excludeDBcol", "Exclude database text column from gene target list", TRUE),
      actionButton("updateTargets", "Update miRNA target list"),
      downloadButton('downloadTargets', 'Download miRNA target list')
    )
  ),
  mainPanel(    
    h4(verbatimTextOutput("datasetName")),
    tabsetPanel(id="tabs1",                   
      tabPanel("Heatmap", plotOutput("heatmapPlot", height=800)), 
      tabPanel("All Wells", plotOutput("scatterPlot", height=800)),
      tabPanel("Plate Viewer", showOutput("scatterPlotI", "dimple"), plotOutput("normalizationComparison", height=800)),
      tabPanel("Hits Plot", showOutput("scatterPlotHits", "dimple")),
      tabPanel("Hits List", chartOutput("table", "datatables")),
      tabPanel("Consensus Hits Plot", showOutput("normcomparison", "polycharts")),
      tabPanel("Consensus Hits List", chartOutput("consensusHitList", "datatables")),
      tabPanel("Consensus Venn Diagram", plotOutput("consensusVennDiagram", height=400)),
      tabPanel("Controls", plotOutput("controlPlot", height=800), plotOutput("rowAndColumn", height=800)),
      tabPanel("miRcancer DB", chartOutput("mircancerTable", "datatables")),
      tabPanel("Target Genes", chartOutput("targets", "datatables")),
      tabPanel("Interaction Table", chartOutput("interactionTable", "datatables")),
      tabPanel("Interaction Graph", plotOutput("interactionGraph")),
      tabPanel("GO Enrichment", chartOutput("goEnrichmentTable", "datatables"))
    )
  )
))

