output$uiOutput_hits_options <- renderUI({
  sidebarPanel(           
    selectInput("normalization", "Raw Data / Normalization:", 
                choices = c("Raw", "poc", "npi", "centered", "rcentered", "zscore", "rzscore", "Bscore", "posEffectNorm")),
    selectInput("method", "Method for margin calculation:", choices = c("SD", "MAD", "quartile")),
    sliderInput("margin", "Margin:",  min = 0.5, max = 5, value = 2.5, step= 0.5),
    checkboxInput("showFilterOptions", "Sample filter options", FALSE),               
    conditionalPanel(
      condition = "input.showFilterOptions",
      helpText("Filter options:"),
      textInput("exclude", "Exclude miRNAs by regular expression:", value=""),
      actionButton("updateExclusion", "Update exclusion filter"),
      textInput("include", "Include miRNAs in hit list by regular expression:", value=""),
      actionButton("updateInclusion", "Update inclusion filter"),
      helpText("For example, you can select all let-7 like this: let-7, or you can select several miRs like this: mir-(765|558)")
    ),
    actionButton("updateNormalization", "Update Settings"))
})

output$uiOutput_hits <- renderUI({
 elements <- list(                         
  tabPanel("Heatmap", fluidRow(               
    column(4,
           checkboxInput("showLabelsOnHeatmap", "Show sample labels in heatmap", TRUE)
    ),
    column(4,
           selectInput("colorA", "Select low signal color for heatmap:", 
                       choices = c("red", "blue", "darkblue", "steelblue", "magenta", "yellow", "white", "green"))
    ),
    column(4, 
           selectInput("colorB", "Select high signal color for heatmap:", c("yellow", "red", "darkblue", "blue", "steelblue", "magenta", "white", "green"))
    )
  ),
  plotOutput("heatmapPlot", height=800)
  ),                
  tabPanel("All Wells", plotOutput("scatterPlot", height=800)),
  tabPanel("Plate Viewer",     
           checkboxInput("showSEM", "Show SEM in bubble charts", TRUE),
           sliderInput("plateSelected", "Select a plate for the plate view:", min=1, max=12, value=1, step=1),
           showOutput("scatterPlotI", "dimple"), 
           plotOutput("normalizationComparison", height=800)
  ),
  tabPanel("Hits Plot", showOutput("scatterPlotHits", "dimple")),
  tabPanel("Hits List", 
           dataTableOutput("table"),
           downloadButton('downloadHits', 'Download hit list')
  ))  
 
 do.call(tabsetPanel, elements)
})