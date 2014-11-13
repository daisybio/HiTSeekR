library(shiny)
require(rCharts)
require(ggplot2)
require(grid)
require(gridExtra)
require(scales)
require(RmiR)

shinyUI(navbarPage(
  HTML('<img src="RNAice.png"/>'),
  windowTitle="RNAice",
  tabPanel("Data",
           conditionalPanel(condition="input.file==null",
             selectInput("dataset", "Select a demo dataset or...", choices = c("BCSC/MaSC miRNA inhibitors" = "BCSC", "MTS data" = "MTS data"))
           ),
    fileInput("file", "Upload a new data set", multiple=FALSE),
    selectInput("fileSeparator", "Column separator", c("tab"= "\t", "comma"= ",", "semicolon"=";")),                     
    uiOutput("uiOutput_data_options"),
    uiOutput("uiOutput_data")
  ),           
  tabPanel("Hits", uiOutput("uiOutput_hits_options"), mainPanel(uiOutput("uiOutput_hits"))),        
  tabPanel("Consensus Hits", uiOutput("uiOutput_consensus_hits")),
  tabPanel("miRNA target genes",uiOutput("uiOutput_mirna_targets")),
  tabPanel("KPM", verbatimTextOutput("KPM.test"), actionButton("startKPMButton", "Start miRNA target gene enrichment with KPM"), downloadButton('downloadIndicatorMatrix')),
  tabPanel("Controls", plotOutput("controlPlot", height=800), plotOutput("rowAndColumn", height=800)),
  tabPanel("miRcancer DB", uiOutput("uiOutput_mircancer")),
  tabPanel("Gene Ontology", uiOutput("uiOutput_geneOntology"))
))