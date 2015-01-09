library(shiny)
require(rCharts)
require(ggplot2)
require(grid)
require(gridExtra)
require(scales)
require(RmiR)
require(shinysky)

shinyUI(fluidPage(
 navbarPage(
  HTML('<img src="RNAice.png"/>'),
  id="mainNavbar",  
  tabPanel("Data",wellPanel(fluidRow(column(4,
           conditionalPanel(condition="input.file==null",
             selectInput("dataset", "Select a demo dataset or...", choices = c("BCSC/MaSC miRNA inhibitors" = "BCSC", "MTS data" = "MTS data"))
           )), column(4,
    fileInput("file", "Upload a new data set", multiple=FALSE)), column(4,
    selectInput("fileSeparator", "Column separator", c("tab"= "\t", "comma"= ",", "semicolon"=";"))))),                     
    checkboxInput("showColOptions", "Show file input options", FALSE),
    hr(),
    conditionalPanel(condition="input.showColOptions", wellPanel(uiOutput("uiOutput_data_options"))),
    uiOutput("uiOutput_data")
  ),   
  tabPanel("Hits", shinyalert("hits_error"), uiOutput("uiOutput_hits_options"), uiOutput("uiOutput_hits")),        
  tabPanel("Consensus Hits", uiOutput("uiOutput_consensus_hits")),
  tabPanel("miRNA target genes",uiOutput("uiOutput_mirna_targets")),
  tabPanel("KPM", verbatimTextOutput("KPM.test"), actionButton("startKPMButton", "Start miRNA target gene enrichment with KPM"), downloadButton('downloadIndicatorMatrix')),
  tabPanel("Controls", plotOutput("controlPlot", height=800), plotOutput("rowAndColumn", height=800)),
  tabPanel("miRcancer DB", uiOutput("uiOutput_mircancer")),
  tabPanel("Gene Ontology", uiOutput("uiOutput_geneOntology"))
)))