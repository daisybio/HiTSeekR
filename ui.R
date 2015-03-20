library(shiny)
require(rCharts)
require(ggplot2)
require(grid)
require(gridExtra)
require(scales)
require(RmiR)
require(shinysky)

shinyUI(navbarPage(
  title=HTML('<img style="height:40px; margin-top: -7px;" src="RNAice.png"/>'),
  windowTitle="RNAice - RNAi comprehensive evaluation",
  header=tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "rnaice.css"),
    shinyalert("general_status")
  ),
  id="mainNavbar",  
  tabPanel("Data",wellPanel(
    tags$style(type="text/css", '#loadOptionsPanel { max-width:800px;}'),
    id="loadOptionsPanel",
    fluidRow(
      column(4,
        conditionalPanel(condition="input.file==null",
             selectInput("dataset", "Select a demo dataset or...", choices = demo.data.sets)
        )), 
      column(4,
            fileInput("file", "Upload a new data set", multiple=FALSE)
        ), 
      column(4,
        selectInput("fileSeparator", "Column separator", c("tab"= "\t", "comma"= ",", "semicolon"=";"))))
        ),                     
    checkboxInput("showColOptions", "Show file input options", FALSE),
    actionButton("startButton", "Process raw data", styleclass="primary"),    
    hr(),
    conditionalPanel(condition="input.showColOptions", 
          wellPanel(
            tags$style(type="text/css", '#dataOptionsPanel { max-width:1200px;}'),
            id="dataOptionsPanel",
            uiOutput("uiOutput_data_options")
          )
    ),
    uiOutput("uiOutput_data")
  ),
  tabPanel("Quality Control", uiOutput("uiOutput_quality_control")),  
  tabPanel("Hit Discovery", shinyalert("hits_error"), uiOutput("uiOutput_hits_options"), uiOutput("uiOutput_hits")),        
  tabPanel("Consensus Hits", uiOutput("uiOutput_consensus_hits")),
  tabPanel("microRNAs", 
           shinyalert("mirna_target_status"), 
           uiOutput("uiOutput_mirna_targets")
  ),
  tabPanel("Systems Biology", uiOutput("uiOutput_geneOntology"))
))