library(shiny)
library(rCharts)
library(ggplot2)
library(grid)
library(gridExtra)
library(scales)
library(RmiR)
library(shinysky)
library(shinyjs)

elts <- list(
  #title=HTML('<img style="height:40px; margin-top: -7px;" src="RNAice.png"/>'),
  #windowTitle="RNAice - RNAi comprehensive evaluation",
  title=HTML('<img style="height:40px; margin-top: -7px;" src="HiTSeekR.png" onclick="window.location.reload(true); 
"/>'),
  windowTitle="High-throughput Screening Kit for R",
  header=tags$head(
    tags$script(
                HTML('
                     Shiny.addCustomMessageHandler("enableNavTab", function(tab){
                      if($(".nav li:nth-child(" + tab + ")").hasClass("disabled")){
                        $(".nav li:nth-child(" + tab + ")").removeClass("disabled");      
                        $(".nav li:nth-child(" + tab + ")").toggle();   
                        $(".nav li:nth-child(" + tab + ")").pulse({times: 4, duration: 1000});      
                      }
                     });
                     Shiny.addCustomMessageHandler("disableNavTab", function(tab){
                      $(".nav li:nth-child(" + tab + ")").addClass("disabled");      
                      $(".nav li:nth-child(" + tab + ")").toggle();
                     });
                     $(document).ready(function(){
                      $(".nav li").attr("class", "disabled");
                      $(".nav li").toggle();
                      $(".nav li").click(function(e){
                        if($(this).hasClass("disabled")){ return false;}
                      });
                      $(".nav li:nth-child(1)").removeClass("disabled");
                      $(".nav li:nth-child(1)").toggle();
                      $(".nav li:nth-child(1)").attr("class", "active");
                      $(".nav li:nth-child(8)").removeClass("disabled");
                      $(".nav li:nth-child(8)").toggle();                      

                     });

    $.fn.pulse = function(options) {
    
        var options = $.extend({
            times: 3,
            duration: 1000
        }, options);
    
        var period = function(callback) {
            $(this).animate({opacity: 0}, options.duration, function() {
                $(this).animate({opacity: 1}, options.duration, callback);
            });
        };
        return this.each(function() {
            var i = +options.times, self = this,
            repeat = function() { --i && period.call(self, repeat) };
            period.call(this, repeat);
        });
    };

    function openOverlay(olEl) {
        $oLay = $(olEl);
                           
                     $("#overlay-shade").fadeTo(1, 1, function() {
                       var props = {
                         oLayWidth       : $oLay.width(),
                         scrTop          : $(window).scrollTop(),
                         viewPortWidth   : $(window).width()
                       };
                       
                       var leftPos = (props.viewPortWidth - props.oLayWidth) / 2;
                       
                       $oLay
                       .css({
                         display : "block",
                         opacity : 0,
                         top : "-=300",
                         left : leftPos+"px"
                       })
                       .animate({
                         top : props.scrTop + 40,
                         opacity : 1
                       }, 600);
                     });
                     }

    function closeOverlay() {
      $(".overlay").animate({
        top : "-=300",
        opacity : 0
      }, 400, function() {
        $("#overlay-shade").fadeOut(300);
        $(this).css("display","none");
      });
    }
                         ')),
    tags$link(rel = "stylesheet", type = "text/css", href = "HiTSeekR.css"),
    shinyalert("general_status")
  ),
  id="mainNavbar",  
  position="fixed-top",  
  tabPanel("Input",
     conditionalPanel("input.dataset == 'none selected'",
       HTML('<div class="shinyalert alert fade alert-info in">Before you can start your analysis you have to select or upload a dataset:</div>'), 
        fluidRow(
          tags$style(type="text/css", '#loadOptionsPanel { max-width:800px;}'),
          id="loadOptionsPanel",
          column(4, wellPanel( tags$style(type="text/css", "#demoPanel { height: 220px; }"),
                               id = "demoPanel",
                 selectInput("screenType", "Type of screen", c("Gene silencing" = "siRNA", "miRNA (e.g. inhibitor)" = "miRNA", "Small compounds (e.g. drugs)" = "compound")),
                 conditionalPanel(condition="input.file==null",
                                  selectInput("dataset", "Select a demo dataset", choices = c("none selected"), "none selected")
                 ))), 
          column(4, wellPanel( tags$style(type="text/css", "#fileLoadPanel { height: 220px; }"),
                               id = "fileLoadPanel",
                 fileInput("file", "Upload a new data set", multiple=FALSE)
          )), 
          column(4, wellPanel( tags$style(type="text/css", "#pubchemloadPanel { height: 220px; }"),
                               id = "pubchemloadPanel",
                 numericInput("pubchem_aid", label = "Use PubChem Assay ID (AID)",value = 743456),     
                 actionButton("getFromPubChemButton", "Download", styleclass="primary")      
          ))
        )),   
    conditionalPanel("input.dataset != 'none selected'",
    HTML('<div class="shinyalert alert fade alert-info in">In the next step you can customize how the data set is processed. Note:
         <ul>
          <li>In particular for uploaded data sets you need to select which type of information is found in which column. Check "show fileinput options".</li>
          <li>log2 transformation is typically recommended for signal data</li>
          <li>The B-score normalization is ideally suited to address position bias in the data, but is computationally expensive and not selected by default.</li>
          <li>If you would like to analyze a different data set just click on the HiTSeekR logo in the top left corner to get back to the start.</li>
          <li>On the bottom of the page, a preview of the input data is shown. Once you are confident about the selected settings press the "Process raw data" button below.</li>
         </ul></div>'), 
    checkboxInput("log2normalize", "Log2 transform signal data", TRUE),
    checkboxInput("computeBscore", "Compute B-score", FALSE),
    conditionalPanel("input.file",
                     selectInput("fileSeparator", 
                                 "Column separator", 
                                 c("tab"= "\t", "comma"= ",", "semicolon"=";"))
    ),
    checkboxInput("showColOptions", "Show file input options", FALSE),
    fluidRow(column(1, actionButton("startButton", "Process raw data", styleclass="primary")),    
    column(1, conditionalPanel("input.startButton",
      actionButton("continueToQC", "Continue with quality control", styleclass="info")
    ))),
    hr(),
    conditionalPanel(condition="input.showColOptions", 
                     wellPanel(
                       tags$style(type="text/css", '#dataOptionsPanel { max-width:1200px;}'),
                       id="dataOptionsPanel",
                       shinyalert("data_processing_status", click.hide = FALSE), 
                       uiOutput("uiOutput_data_options")
                     )
    ),
    tabsetPanel(
      tabPanel("Preview of input data",
        dataTableOutput("table_rawData")
      )
    ))
  ),  
  tabPanel("Quality Control", id="qc", 
           HTML('<div class="shinyalert alert fade alert-info in">Below you can select plots that may help to assess potential quality problems with the raw data.</div>'), 
           actionButton("continueToNormalizationEffect", "Continue with studying normalization effects", styleclass="info"), checkboxInput("showHelpPages", "Show help text", FALSE), uiOutput("uiOutput_quality_control")),
  tabPanel("Normalization Effect", 
           HTML('<div class="shinyalert alert fade alert-info in">Some quality issues can be accommodated with appropriate normalization. Here you can study the effect of different normalization strategies on the data:</div>'), 
           actionButton("continueToHitDiscovery", "Continue with hit discovery", styleclass="info"),
           hr(),
           fluidRow(column(width=4, uiOutput("uiOutput_dataWellPanel"))),uiOutput("uiOutput_data")),
  tabPanel("Hit Discovery", shinyalert("hits_error"), 
           HTML('<div class="shinyalert alert fade alert-info in">Here you can select a normalization and hit detection strategy to generate a candidate hit list used for down-stream analysis :</div>'), 
           conditionalPanel("input.screenType == 'miRNA'",
                            actionButton("continueToMiRNAs", "Continue with miRNA analysis and target discovery", styleclass="info")           
           ),
           conditionalPanel("input.screenType == 'compound'",
                            actionButton("continueToDrugs", "Continue with drug target discovery", styleclass="info")           
           ),
           conditionalPanel("input.screenType == 'siRNA'",
                            actionButton("continueToGenes", "Continue with gene list analyses", styleclass="info")
           ),
           hr(),
           uiOutput("uiOutput_hits_options"), uiOutput("uiOutput_hits")),        
  #tabPanel("Consensus Hits", 
  #         uiOutput("uiOutput_consensus_hits_options"), 
  #         mainPanel(uiOutput("uiOutput_consensus_hits"))
  #),
  tabPanel("microRNAs",            
           shinyalert("mirna_target_status"), 
           HTML('<div class="shinyalert alert fade alert-info in">Here you can investigate miRNA targets, family membership and find publications associated with the hit candidates found in the previous tab. the effect of different normalization strategies on the data:</div>'), 
           actionButton("continueToGenes2", "Continue with gene list analyses", styleclass="info"),
           hr(),
           uiOutput("uiOutput_mirna_targets")
  ),
  tabPanel("Small Compounds", 
           HTML('<div class="shinyalert alert fade alert-info in">Here you can find putative drug target proteins of the previously identified hit candidates.</div>'), 
           actionButton("continueToGenes3", "Continue with gene list analyses", styleclass="info"),
           hr(),
           uiOutput("uiOutput_drug_targets")
  ),
  tabPanel("Genes", 
           HTML('<div class="shinyalert alert fade alert-info in">Here you can perform down-stream analysis of genes identified in the previous step:</div>'), 
           uiOutput("uiOutput_gene_set_analysis")
  ),
  tabPanel("About",
           fluidRow(column(width=4, HTML("The High-Throughput Screening kit for R (HiTSeekR) was developed as a joint project between<br/><br/>
                <div style='float:left;'><a href='http://nanocan.org'><img width=150 src='NanoCAN.png'/></a></div>
                <div style='float:right;'><a href='http://baumbachlab.net'><img width=200 src='baumbachlab.png'/></a></div>
                <div style='clear: both; padding-top:50px;'>
                Contact: Markus List &lt;mlist'at'health.sdu.dk&gt;
                </div>
                ")))
  )
)

shinyUI(fluidPage(
    HTML('<div id="overlay-shade"></div>
          <div id="overlay-inAbox" class="overlay">
          <div class="wrapper">
            <img style="height:40px; margin: 50px;" src="HiTSeekR.png"/>
            <p>Welcome to the High-Throughput Screening Kit for R. 
This web application is dedicated to the analysis of high-throughput screening data of various types. It can accommodate small to ultra-large scale. Start by selecting a screen type below.</p>
          </div>
          <div class="toolbar">                
         '),
    HTML('If you are here the first time, check out the <a target="_blank" href="http://nanocan.github.io/HiTSeekR/"><button id="tutorial" type="button" class="btn action-button btn-info shiny-bound-input">Tutorial
    </button></a>'),
    h1("Select type of screen:"),
    actionButton("siRNA", "Gene silencing", "primary"),
    actionButton("miRNA", "MicroRNA", "primary"),
    actionButton("compound", "Small compound", "primary"),
    HTML('</div></div>
         <script>    openOverlay("#overlay-inAbox"); </script>'),
    useShinyjs(),
    do.call(navbarPage, elts)    
  )
)

# shinyUI(fluidRow(
#   column(12,
#          "", 
#          do.call(navbarPage, elts)
#   )
# ))