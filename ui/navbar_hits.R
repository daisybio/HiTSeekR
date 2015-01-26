#return experiments included in the selected data set
experiments <- reactive({
  experiments <- unique(as.character(processedData()$Experiment))
  if(is.null(experiments)) return(NULL)
  else return(experiments)
})

output$uiOutput_hits_options <- renderUI({
  wellPanel(
  fluidRow(
    column(4,
      actionButton("updateNormalization", "Update Settings", styleclass="primary"),
      checkboxInput("show.sem.in.hits", "Show replicate standard error of the mean", FALSE)
    ),column(4,      
      selectInput("experimentSelected", "Select experiment:", experiments(), experiments()[1]),
      selectInput("normalization", "Raw Data / Normalization:", 
                choices = c("Raw signal" = "Raw", "Percentage of negative control (poc)" = "poc", "Normalized percentage inhibition (npi)" = "npi", "Centered by mean (centered)" = "centered", "Centered by median (rcentered)" = "rcentered", "z-score" = "zscore", "Robust z-score (rzscore)" = "rzscore", "B-score" = "Bscore"))      
    ),column(4,
      selectInput("method", "Method for margin calculation:", choices = c("Bayesian Statistics" = "Bayes", "Standard deviation" ="SD", "Median absolute deviation" = "MAD", "Inter-quartile range" = "quartile")),
      conditionalPanel(
        condition = "input.method != 'Bayes'",
        sliderInput("margin", "Margin:",  min = 0.5, max = 5, value = 2.5, step= 0.5)
      ),
      conditionalPanel(
        condition = "input.method == 'Bayes'",
        wellPanel(
        selectInput("bayes_hypothesis", "Select hypothesis:", 
                    choices = list("effect" = "p_effect", "suppressor" = "p_suppressor", "promotor"= "p_promotor")),
        sliderInput("bayes_pval_cutoff", "p-value cutoff (>)", min=0, max=1, value = 0.95, step=0.01)
        )
      )
    )
  ), fluidRow(
    checkboxInput("showFilterOptions", "Sample filter options", FALSE),               
    conditionalPanel(
      condition = "input.showFilterOptions",
      helpText("Filter options:"),
      textInput("exclude", "Exclude miRNAs by regular expression:", value=""),
      actionButton("updateExclusion", "Update exclusion filter"),
      textInput("include", "Include miRNAs in hit list by regular expression:", value=""),
      actionButton("updateInclusion", "Update inclusion filter"),
      helpText("For example, you can select all let-7 like this: let-7, or you can select several miRs like this: mir-(765|558)")
    ))
 )
})

output$uiOutput_hits <- renderUI({
  exp.data <- processedData()
  plates <- as.integer(unique(exp.data$Plate))
  replicates <- as.integer(unique(exp.data$Replicate))
  
 elements <- list(    
   #tabPanel("Scatter Plots", uiOutput("scatterPlotTagList")),
   tabPanel("Plate Viewer",
            sidebarPanel(
              sliderInput("plateSelected", "Select a plate:", min=min(plates), max=max(plates), value=min(plates), step=1),
              sliderInput("replicateSelected", "Select a replicate (heatmap):", min=min(replicates), max=max(replicates), value=min(replicates), step=1)
            ),mainPanel(
              tags$script(src = "https://code.highcharts.com/modules/heatmap.js"),
              showOutput("intHeatmapPlot", "highcharts"),
              showOutput("intPlateScatterPlot", "highcharts")
              #showOutput("scatterPlotI", "dimple")
            #plotOutput("normalizationComparison", height=800
            )
   ),   
  tabPanel("Heatmap", fluidRow(               
           checkboxInput("showLabelsOnHeatmap", "Show sample labels in heatmap", FALSE)
  ),
  plotOutput("heatmapPlot", height=800)
  ),                
  tabPanel("Hits Plot", tags$div(showOutput("scatterPlotHits", "dimple"))),
  tabPanel("Hits List", 
           dataTableOutput("table_hits"),
           downloadButton('downloadHits', 'Download hit list')
  ),
  tabPanel("miRNA family coverage", 
           shinyalert("miRNA_family_info", click.hide = TRUE),
           wellPanel(fluidRow(
            column(6, sliderInput("family_size_cutoff", "Family size cutoff:", min=0, max=20, value=0)),
            column(6, sliderInput("family_coverage_cutoff", "Family coverage cutoff:", min=0, max=100, value=0))
           )),
           dataTableOutput("family.hitrate"))
 )  
 
 do.call(tabsetPanel, elements)
})