library(shiny)
library(rCharts)
library(ggplot2)
library(grid)
library(gridExtra)
library(scales)
library(RmiR)
library(shinysky)

elts <- list(
  #title=HTML('<img style="height:40px; margin-top: -7px;" src="RNAice.png"/>'),
  #windowTitle="RNAice - RNAi comprehensive evaluation",
  title=HTML('<img style="height:40px; margin-top: -7px;" src="HiTSeekR.png"/>'),
  windowTitle="High-throughput Screening Kit for R",
  header=tags$head(
    tags$script(
                HTML('
                     Shiny.addCustomMessageHandler("enableNavTab", function(tab){
                      $(".nav li:nth-child(" + tab + ")").removeClass("disabled");      
                     });
                     Shiny.addCustomMessageHandler("disableNavTab", function(tab){
                      $(".nav li:nth-child(" + tab + ")").addClass("disabled");      
                     });
                     $(document).ready(function(){
                      $(".nav li").attr("class", "disabled");
                      $(".nav li").click(function(e){
                        if($(this).hasClass("disabled")){ return false;}
                      });
                      $(".nav li:nth-child(1)").removeClass("disabled");
                      $(".nav li:nth-child(1)").attr("class", "active");
                     });
                     ')),
    #tags$link(rel = "stylesheet", type = "text/css", href = "rnaice.css"),
    #tags$link(rel = "stylesheet", type = "text/css", href = "supernice.css"),
    tags$link(rel = "stylesheet", type = "text/css", href = "HiTSeekR.css"),
    shinyalert("general_status")
  ),
  id="mainNavbar",  
  tabPanel("Input",wellPanel(
    tags$style(type="text/css", '#loadOptionsPanel { max-width:800px;}'),
    id="loadOptionsPanel",
    fluidRow(
      column(4,
             conditionalPanel(condition="input.file==null",
                              selectInput("dataset", "Select a demo dataset", choices = demo.data.sets)
             )), 
      column(4,
             fileInput("file", "Upload a new data set", multiple=FALSE),
             selectInput("fileSeparator", 
                         "Column separator", 
                         c("tab"= "\t", "comma"= ",", "semicolon"=";"))
      ), 
      column(4,
             numericInput("pubchem_aid", label = "PubChem AID",value = 743456),     
             actionButton("getFromPubChemButton", "Download from PubChem")      
      )
    )),                     
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
    dataTableOutput("table_rawData")
  ),  
  tabPanel("Data", uiOutput("uiOutput_dataWellPanel"),uiOutput("uiOutput_data")),
  tabPanel("Quality Control", id="qc", uiOutput("uiOutput_quality_control")),
  tabPanel("Hit Discovery", shinyalert("hits_error"), uiOutput("uiOutput_hits_options"), uiOutput("uiOutput_hits")),        
  tabPanel("Consensus Hits", 
           uiOutput("uiOutput_consensus_hits_options"), 
           mainPanel(uiOutput("uiOutput_consensus_hits"))
  ),
  tabPanel("microRNAs", 
           shinyalert("mirna_target_status"), 
           uiOutput("uiOutput_mirna_targets")
  ),
  tabPanel("Small Compounds", 
           uiOutput("uiOutput_drug_targets")
  ),
  tabPanel("Genes", 
           uiOutput("uiOutput_gene_set_analysis")
  )
)


shinyUI(fluidRow(
  column(12,
         "", 
         do.call(navbarPage, elts)
  )
))