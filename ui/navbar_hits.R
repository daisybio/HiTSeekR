output$uiOutput_hits_options <- renderUI({
 wellPanel(
  fluidRow(
    column(4,
      actionButton("updateNormalization", "Update Settings", styleclass="primary")
    ),column(4,      
      selectInput("experimentSelected", "Select experiment:", unique(as.character(processedData()$Experiment)), processedData()$Experiment[1]),
      selectInput("normalization", "Raw Data / Normalization:", 
                choices = c("Raw", "poc", "npi", "centered", "rcentered", "zscore", "rzscore", "Bscore", "posEffectNorm"))
    ),column(4,
      selectInput("method", "Method for margin calculation:", choices = c("SD", "MAD", "quartile")),
      sliderInput("margin", "Margin:",  min = 0.5, max = 5, value = 2.5, step= 0.5)
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


output$scatterPlotTagList <- renderUI({
  exp.data <- processedData()
  plots <- foreach(experiment = unique(exp.data$Experiment)) %do% {
    plot.name <- paste(experiment, "IntScatterPlot", sep="")
    showOutput(plot.name, "highcharts")
  }
  do.call(tagList, plots)
})

output$uiOutput_hits <- renderUI({
  exp.data <- processedData()
  plates <- as.integer(unique(exp.data$Plate))
  replicates <- as.integer(unique(exp.data$Replicate))
  
 elements <- list(    
   tabPanel("Scatter Plots", uiOutput("scatterPlotTagList")),
   tabPanel("Plate Viewer",
            sidebarPanel(
              checkboxInput("showSEM", "Show SEM", TRUE),
              sliderInput("plateSelected", "Select a plate:", min=min(plates), max=max(plates), value=min(plates), step=1),
              sliderInput("replicateSelected", "Select a replicate:", min=min(replicates), max=max(replicates), value=min(replicates), step=1)
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