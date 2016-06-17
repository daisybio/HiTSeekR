#return experiments included in the selected data set
experiments <- reactive({
  experiments <- unique(as.character(processedData()$Experiment))
  if(is.null(experiments)) return(NULL)
  else return(experiments)
})

#return readouts included in the selected data set
readouts <- reactive({
  readouts <- unique(as.character(processedData()$Readout))
  if(is.null(readouts)) return(NULL)
  else return(readouts)
})

minCutoff <- reactive({
  exp.data <- processedData()
  
  if(is.null(exp.data)) return(NULL)
  return(round(min(exp.data[["Raw"]], na.rm=T)))
})

maxCutoff <- reactive({
  exp.data <- processedData()
  
  if(is.null(exp.data)) return(NULL)
  return(round(max(exp.data[["Raw"]], na.rm=T)))
})

defaultMaxCutoff <- reactive({
  minCut <- minCutoff()
  maxCut <- maxCutoff()
  
  if(is.null(minCut) || is.null(maxCut)) return(NULL)
  else return(round((maxCut - minCut) * 0.8 + minCut))
})

defaultMinCutoff <- reactive({
  minCut <- minCutoff()
  maxCut <- maxCutoff()
  
  if(is.null(minCut) || is.null(maxCut)) return(NULL)
  else return(round((maxCut - minCut) * 0.2 + minCut))
})

normalizationChoices <- reactive({
  methods <- c("Raw signal" = "Raw", "Centered by mean (centered)" = "centered", "Centered by median (rcentered)" = "rcentered", "z-score" = "zscore", "Robust z-score (rzscore)" = "rzscore")
  controlBasedMethods <- c("Percentage of negative control (poc)" = "poc", "Normalized percentage inhibition (npi)" = "npi")
  Bscore <- c("B-score" = "Bscore")
  if(input$computeBscore){
    methods <- c(methods, Bscore)
  }
  if(input$hasControls){
    if(!is.null(input$posCtrl))
      return(c(methods, controlBasedMethods))
    else{
      return(c(methods, controlBasedMethods[1])) 
    }
  }
  else return(methods)
})

marginChoices <- reactive({
  methods <- c("Standard deviation" ="SD", "Median absolute deviation" = "MAD", "Inter-quartile range" = "quartile", "Fixed cutoffs" = "cutoff")
  controlBasedMethods <- c("Bayesian Statistics" = "Bayes", "SSMD" = "SSMD")
  if(input$hasControls)
    return(c(methods, controlBasedMethods))
  else return(methods)
})

output$uiOutput_hits_options <- renderUI({
  wellPanel(
    tags$style(type="text/css", '#hitsOptionsPanel { max-width:1400px;}'),
    id="hitsOptionsPanel",
    fluidRow(
      column(3,
             checkboxInput("showSamplePosition", "Show sample position", FALSE),      
             checkboxInput("showAllScores", "Show all computed scores", FALSE),      
             conditionalPanel(condition = "input.replicateCol && input.showAllScores",
                              checkboxInput("show.sd.in.hits", "Show replicate standard deviation", FALSE)
             ),
             checkboxInput("showFilterOptions", "Sample filter options", FALSE),
             checkboxInput("differentialScreening", "Differential screening", FALSE)
      ),column(2,
               selectInput("experimentSelected", "Select experiment:", experiments(), experiments()[1], multiple=TRUE),
               selectInput("readoutSelected", "Select readout:", readouts(), readouts()[1], multiple=TRUE)
      ),column(2,
               selectInput("normalization", "Normalization:", 
                           choices = normalizationChoices()),      
               selectInput("method", "Hit detection method:", choices = marginChoices())      
      ),column(2,
               selectInput("effect", "Effect:", choices = c("effect", "suppressor", "promotor")),
               conditionalPanel(
                 condition = "input.method != 'Bayes' && input.method != 'cutoff'",
                 sliderInput("margin", "Margin:",  min = 0, max = 20, value = 3, step= 0.5)
               ),
               conditionalPanel(
                 condition = "input.method == 'cutoff'",
                 numericInput("upperCutoff", "Upper cutoff:", value=defaultMaxCutoff()),
                 numericInput("lowerCutoff", "Lower cutoff:", value=defaultMinCutoff())
               ),
               conditionalPanel(
                 condition = "input.method == 'Bayes'",
                 wellPanel(
                   sliderInput("bayes_pval_cutoff", "p-value cutoff (>)", min=0, max=1, value = 0.95, step=0.01),
                   actionButton("computeBayes", "Apply Bayes method", styleclass="primary")
                 )
               )
      ),column(3, uiOutput("cat_legend"))
    ),
    fluidRow(
      column(4,
             conditionalPanel(
               condition = "input.differentialScreening",
               helpText("Differential screening options:"),
               conditionalPanel(
                 condition = "input.method != 'Bayes'",
                 sliderInput("diffMargin", "Reference margin:",  min = 0.5, max = 10, value = 2.5, step= 0.5)
               ),
               selectInput("referenceExperiment", "Reference experiment:", experiments(), experiments()[1]),
               selectInput("referenceReadout", "Reference readout:", readouts(), readouts()[1]),         
               helpText("Hits found in the reference readout (experiment) are removed from the target readouts (experiments)")
             ) 
      )),
    fluidRow(
      column(4,
             conditionalPanel(
               condition = "input.showFilterOptions",
               helpText("Filter options:"),
               textInput("exclude", "Exclude sample in hit list by regular expression:", value=""),
               textInput("include", "Include sample in hit list by regular expression:", value=""),    
               helpText("For example, you can select all let-7 like this: let-7, or you can select several miRs like this: mir-(765|558)")
             ) 
      ))
  )
})

output$uiOutput_hits <- renderUI({
  exp.data <- processedData()
  plates <- as.integer(unique(exp.data$Plate))
  replicates <- as.integer(unique(exp.data$Replicate))
  
  elements <- list(    
    #tabPanel("Scatter Plots", uiOutput("scatterPlotTagList")),
    tabPanel("Hits List", 
             conditionalPanel("input.showHelpText",
                              HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                         hitsTableInfoText,
                                         '</div>', sep="")
                              )
             ),
             dataTableOutput("table_hits"),
             downloadButton('downloadHits', 'Download hit list')
    ),
    tabPanel("Hits Plot", 
             conditionalPanel("input.showHelpText",
                              HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                         hitsPlotInfoText,
                                         '</div>', sep="")
                              )
             ),showOutput("scatterPlotHits", "dimple")),
    tabPanel("Heatmap", 
             conditionalPanel("input.showHelpText",
                              HTML(paste('<div class="shinyalert alert fade alert-info in">',
                                         heatmapInfoText,
                                         '</div>', sep="")
                              )
             ),
             fluidRow(  
               selectInput("heatmapExperimentSelected", "Select experiment:", input$experimentSelected),
               selectInput("heatmapReadoutSelected", "Select readout:", input$readoutSelected),         
               checkboxInput("showLabelsOnHeatmap", "Show sample labels in heatmap", FALSE)
             ),
             plotOutput("heatmapPlot", height=800)
    ))                 
  do.call(tabsetPanel, elements)
})