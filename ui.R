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
  titlePanel(HTML('<img src="RNAice.png"/>'), windowTitle="RNAice - RNAi comprehensive evaluation"),
  id="mainNavbar",  
  tabPanel("Data",wellPanel(fluidRow(column(4,
           conditionalPanel(condition="input.file==null",
             selectInput("dataset", "Select a demo dataset or...", choices = c("BCSC/MaSC miRNA inhibitors" = "BCSC", "Melanoma-Inhibiting miRNAs" = "A375_MTS"))
           )), column(4,
    fileInput("file", "Upload a new data set", multiple=FALSE)), column(4,
    selectInput("fileSeparator", "Column separator", c("tab"= "\t", "comma"= ",", "semicolon"=";"))))),                     
    checkboxInput("showColOptions", "Show file input options", FALSE),
    hr(),
    conditionalPanel(condition="input.showColOptions", wellPanel(uiOutput("uiOutput_data_options"))),
    uiOutput("uiOutput_data")
  ),   
  tabPanel("Hit Discovery", shinyalert("hits_error"), uiOutput("uiOutput_hits_options"), uiOutput("uiOutput_hits")),        
  tabPanel("Consensus Hits", uiOutput("uiOutput_consensus_hits")),
  tabPanel("miRNA Targets", 
           shinyalert("mirna_target_status"), 
           uiOutput("uiOutput_mirna_targets")),
  tabPanel("Controls", plotOutput("controlPlot", height=800), plotOutput("rowAndColumn", height=800)),  
  tabPanel("Gene Ontology", uiOutput("uiOutput_geneOntology"))
)))