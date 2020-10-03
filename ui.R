library(shiny)
library(rCharts)
library(ggplot2)
library(grid)
library(gridExtra)
library(scales)
library(shinysky)

elts <- list(
  title=HTML('<img style="height:40px; margin-top: -7px;" src="HiTSeekR.png" onclick="window.location.reload(true); 
"/>'),
  windowTitle="HiTSeekR",
  header=tags$head(
    tags$script(
                HTML('
                     Shiny.addCustomMessageHandler("disableScreenType", function(empty){
                        $("#screenType").attr({
                            "disabled": "disabled"
                        });
                        closeOverlay();
                     });
                     Shiny.addCustomMessageHandler("disableLoadingScreen", function(empty){
                        $("#loadingBar").hide();
                        $("#buttonBar").show();
                     });
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
    shinyalert("general_status"), fluidRow(id='helpCheckBox', tags$style(type="text/css", "#helpCheckBox { margin-left:50px;}"), checkboxInput("showHelpText", "Show help", FALSE)),
    HTML('<link rel="stylesheet" type="text/css" href="cookieconsent.min.css"/><script src="cookieconsent.min.js"></script><script>window.addEventListener("load", function(){window.wpcc.init({"colors":{"popup":{"background":"#cff5ff","text":"#000000","border":"#5e99c2"},"button":{"background":"#5e99c2","text":"#ffffff"}}, "padding":"none","margin":"none","fontsize":"tiny","content":{"href":"https://www.learn-about-cookies.com/"},"position":"bottom-right"})});</script>')
  ),
  id="mainNavbar",  
  position="fixed-top",  
  tabPanel("Input",
     conditionalPanel("input.dataset == 'none selected'",
       conditionalPanel("input.showHelpText",                      
        HTML('<div class="shinyalert alert fade alert-info in">Before you can start your analysis you have to upload a dataset (in csv or excel format) or select a demo dataset.</div>')
       ),
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
    conditionalPanel("input.showHelpText",                     
    HTML('<div class="shinyalert alert fade alert-info in">In the next step you can customize how the data set is processed. If questions arise check <a href="https://biomedical-big-data.de/HiTSeekR/tutorial/">tutorial</a>. Note:
         <ul>
          <li>In particular for uploaded data sets you need to select which type of information is found in which column. Check "show fileinput options".</li>
          <li>log2 transformation is typically recommended for signal data</li>
          <li>The B-score normalization is ideally suited to address position bias in the data, but is computationally expensive and not selected by default.</li>
          <li>If you would like to analyze a different data set just click on the HiTSeekR logo in the top left corner to get back to the start.</li>
          <li>On the bottom of the page, a preview of the input data is shown. Once you are confident about the selected settings press the "Process raw data" button below.</li>
         </ul></div>')
    ),
    shinyalert("demoDatasetInfo", click.hide = FALSE),
    conditionalPanel("output.fileUploaded",
                     checkboxInput("isHitList", "If the uploaded file is already a processed list of hits subject to down-stream analysis check this box", FALSE),
                     selectInput("fileSeparator", 
                                 "Column separator", 
                                 c("tab"= "\t", "comma"= ",", "semicolon"=";"))
    ),
    conditionalPanel("!input.isHitList",
                     checkboxInput("log2normalize", "Log2 transform signal data", TRUE),
                     checkboxInput("computeBscore", "Compute B-score", FALSE)
    ),
    checkboxInput("showColOptions", "Show file input options", FALSE),
    conditionalPanel("!input.isHitList", fluidRow(
      tags$style(type="text/css", '#loadButtonsPanel { max-width:400px;}'),
      id="loadButtonsPanel",
      column(6, actionButton("startButton", "Process raw data", styleclass="primary")),    
      column(6, conditionalPanel("input.startButton",
        actionButton("continueToQC", "Continue with quality control", styleclass="info")
    )))),
    conditionalPanel("input.isHitList", fluidRow(
      column(6, actionButton("analysisButton", "Continue with downstream analysis", styleclass="info"))    
    )),
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
           conditionalPanel("input.showHelpText",
             HTML('<div class="shinyalert alert fade alert-info in">Below you can select plots that may help to assess potential quality problems with the raw data.</div>')
           ),
            actionButton("continueToNormalizationEffect", "Continue with studying normalization effects", styleclass="info"),
            hr(),
            uiOutput("uiOutput_quality_control")),
  tabPanel("Normalization Effect", 
           conditionalPanel("input.showHelpText",
            HTML('<div class="shinyalert alert fade alert-info in">Some quality issues can be accommodated with appropriate normalization. Here you can study the effect of different normalization strategies on the data:</div>')
           ),
           fluidRow(
             column(width=4, uiOutput("uiOutput_dataWellPanel")),
             column(4, actionButton("continueToHitDiscovery", "Continue with hit discovery", styleclass="info"))
           ), uiOutput("uiOutput_data")),
  tabPanel("Hit Discovery", 
           conditionalPanel("input.showHelpText",
            HTML('<div class="shinyalert alert fade alert-info in">Here you can select a normalization and hit detection strategy to generate a candidate hit list used for down-stream analysis :</div>')
           ),
           shinyalert("hits_error"),
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
           conditionalPanel("input.showHelpText",
            HTML('<div class="shinyalert alert fade alert-info in">Here you can investigate miRNA targets, miRNA family membership and browse the mircancer database.</div>')
           ),
           actionButton("continueToGenes2", "Continue with gene list analyses", styleclass="info"),
           hr(),
           uiOutput("uiOutput_mirna_targets")
  ),
  tabPanel("Small Compounds", 
           conditionalPanel("input.showHelpText",
             HTML('<div class="shinyalert alert fade alert-info in">Here you can find putative drug target proteins of the previously identified hit candidates.</div>')
           ),
           actionButton("continueToGenes3", "Continue with gene list analyses", styleclass="info"),
           hr(),
           uiOutput("uiOutput_drug_targets")
  ),
  tabPanel("Genes",
           conditionalPanel("input.showHelpText",
            HTML('<div class="shinyalert alert fade alert-info in">Here you can perform down-stream analysis of genes identified in the previous step:</div>')
           ),
           uiOutput("uiOutput_gene_set_analysis")
  ),
  tabPanel("About",
           fluidRow(column(width=8, HTML("The High-Throughput Screening kit for R (HiTSeekR) was developed as a joint project between<br/><br/>
                <a href='http://nanocan.org'><img width=150 src='NanoCAN.png'/></a>
                <a href='https://exbio.de'><img width=200 src='baumbachlab.png'/></a>
                <a href='https://biomedical-big-data.de'><img width=300 src='biomedbigdata_logo.png'/></a>
                <div style='clear: both; padding-top:50px;'>
                Contact: Markus List &lt;markus.list=.AT.=wzw.tum.de&gt;<br/><br/>
                <p><b>If you find HiTSeekR useful in your research please cite:</b></p>
                <p>List, M., Schmidt, S., Christiansen, H., Rehmsmeier, M., Tan, Q., Mollenhauer, J., & Baumbach, J. (2016). Comprehensive analysis of high-throughput screens with HiTSeekR. Nucleic acids research, 44(14), 6639-6648.</a></p>
                <p><a href='https://doi.org/10.1093/nar/gkw554' target='_blank' class='btn btn-primary'>Read the paper &rarr;</a></p>
                <p><a href='https://github.com/biomedbigdata/HiTSeekR' target='_blank' class='btn btn-primary'>View on GitHub &rarr;</a></p>
                <p><a href='https://biomedical-big-data.de/HiTSeekR/' target='_blank' class='btn btn-primary'>Project page with tutorial &rarr;</a></p>
                <br/><br/><hr>
                <div id='aboutdbs'><h4>The following is a list of R packages used in HiTSeekR for annotation and systems biology analysis:</h4>
                <ul>
                <li> reactome.db v. 1.54.1 </li>
                <li> KEGG.db v. 3.2.2 </li>
                <li> GO.db v. 3.2.2</li>
                <li> HTSanalyzeR v. 2.22.0</li>
                <li> mirbase.db v. 1.2.0</li>
                <li> RmiR v. 1.26.0</li>
                </ul>
                <h4>In addition, HiTSeekR integrates the following external resources:</h4>
                <ul>
                  <li><a href='http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=site/page&view=software' target='_blank' class='btn btn-primary'>DIANA tools (web service)</a> &rarr; microRNA target and pathway prediction</li>
                  <li><a href='http://stitch.embl.de/'target='_blank' class='btn btn-primary'>STITCH (v.5.0)</a> &rarr; interaction between proteins and chemicals</li>
                  <li><a href='http://keypathwayminer.compbio.sdu.dk/'target='_blank' class='btn btn-primary'>KeyPathwayMinerWeb</a> &rarr; de novo network enrichment</li>
                </ul>
                </div>
                </div>
                ")))
  )
)

shinyUI(fluidPage(
    HTML('<div id="overlay-shade"></div>
          <div id="overlay-inAbox" class="overlay">
          <div class="wrapper">
            <div style="float:right;">
<!--<img style="height:120px; margin-top: 40px;" src="pipetting.jpg"/>-->
<iframe style="margin-top:50px;" width="180" height="120" src="https://www.youtube.com/embed/CVElY1us2Mw" frameborder="0" allowfullscreen></iframe>
</div>
            <div style="float:left;"><img style="height:40px; margin-left: 10px; margin-bottom: 20px;" src="HiTSeekR.png"/>

<p style="width: 300px; text-align:justify;">Welcome to the High-Throughput Screening Kit for R. 
This web application is dedicated to the analysis of high-throughput screening data of various types. It can accommodate small to ultra-large scale. 
Note that uploaded screening data is only saved transiently during the analysis and will not be disclosed. Start by selecting a screen type below.</p>
          </div></div>
          <div class="toolbar" style="clear:both;">                
         '),
    HTML('<span style="margin-left:20px;"><a target="_blank" href="https://biomedical-big-data.de/HiTSeekR/"><button id="tutorial" type="button" class="btn action-button btn-info shiny-bound-input">If you are here the first time, check out the tutorial
    </button></a></span>'),
    br(),
    HTML("<div id='buttonBar' style='display: none; background-color:#ccccff; padding:30px; margin-top:30px; border-radius: 5px; '><h1>Select type of screen</h1>"),
    br(),
    actionButton("siRNA", "Genes", "primary"),
    actionButton("miRNA", "microRNAs", "primary"),
    actionButton("compound", "Small compounds", "primary"),
    HTML('</div>
         <div id="loadingBar" style="padding:12px; margin-top:30px;">
            <h3>Loading HiTSeekR</h3><p>Please be patient</p> <img src="loading.gif"/>
         </div>
         </div></div>
         <script>   
            openOverlay("#overlay-inAbox"); 
         </script>'),
    do.call(navbarPage, elts)    
  )
)