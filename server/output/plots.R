#dummy progress bar
dummyProgressBar <- function(){
  # Create a Progress object
  progress <- shiny::Progress$new()
  progress$set(message = "Generating plot", value = 0.1)
  # Close the progress when this reactive exits (even if there's an error)
  return(progress)
}

#format integers for plots
formatIntegerForPlot <- function(data, columnToFormat){    
    maxNumber <- max(as.integer(data[,columnToFormat]))
    data[,columnToFormat] <- sprintf(paste("%0", nchar(maxNumber), "d", sep=""),as.integer(data$Plate))
    data[,columnToFormat] <- as.factor(data[,columnToFormat])
    return(data)
}

# bar chart to visualize how many normalizations report a given hit #
output$normcomparison <- renderChart({
  exp.data <- outliers.all()
  consensus <- consensusHitList()
  
  exp.data <- filter(exp.data, Sample %in% consensus$Sample)
  p3 <- rPlot(x="Sample", y="n", color='method', data = exp.data, type = 'bar')
  p3$addParams(dom='normcomparison', width="800", height="600")
  return(p3)
})

# Venn diagram to show overlap of hits in normalization methods #
output$consensusVennDiagram <- renderPlot({
  data <- outliers.all()
  data <- split(data, data$method)
  if(length(data) <= 5){
    data <- lapply(data, function(x){na.omit(x$Sample)})
    fill.colors <- c("purple", "cornflowerblue", "green", "yellow", "orange", "darkorchid1")
    fill.colors <- fill.colors[1:length(data)]
    p <- venn.diagram(data, filename = NULL, force.unique=T, col="transparent", 
                      fill = fill.colors,
                      alpha=0.5, fontfamily="serif")
    grid.draw(p)
  }
  else textplt("this feature only supports up to 5 normalization methods")
})

# Heatmap #
output$heatmapPlot <- renderPlot({  
  plot.data <- data()
  plot.data <- subset(plot.data, Experiment == input$experimentSelected  & Readout == input$readoutSelected)
  hits <- outliers()
  
  my.heatmap(plot.data, input$normalization, input$margin, input$method, hits, showSampleLabels=input$showLabelsOnHeatmap, signalColumn=input$normalization)
})

output$intHeatmapPlot <- renderChart2({
  plot.data <- processedData()
  plot.data <- subset(plot.data, Plate == input$plateSelected & Experiment == input$experimentSelected  & Replicate == input$replicateSelected)
  plot.data <- plot.data[,c("Column", "Row", input$normalization, "Sample", "Accession")]
  colnames(plot.data)[1:3] <- c("x", "y", "value")
  plot.data <- arrange(plot.data, x, y)
  p <- highcharts.heatmap(as.data.frame(plot.data), input$normalization, useWithShiny = T)  
  
  return(p)
})

# Interactive plots for screen hits #
output$scatterPlotHits <- renderChart({
  outl <- outliers()
  p1 <- dPlot(y=input$normalization, x=c("Sample", "Well.position"), z=paste(input$normalization, "_sem", sep=""), 
              data=outl, type="bubble", groups="category", height="800", width="100%",
              bounds = list(x=70, y=30, height="600", width="90%"))  
  p1$addParams(dom='scatterPlotHits')
  #in case we don't see any promotors we need to change the default color to red to make sure suppressors remain colored red
  if(!("promotor" %in% as.character(outl$category)))
    p1$defaultColors("#!d3.scale.ordinal().range(['#FB8072','#FDB462']).domain([0,1])!#")
  else if("included" %in% as.character(outl$category))
    p1$defaultColors("#!d3.scale.ordinal().range(['#FDB462','#80B1D3','#FB8072']).domain([0,1])!#")
  if(input$show.sem.in.hits) p1$zAxis(type="addMeasureAxis", overrideMax = 2*max(outl[[paste(input$normalization, "_sem", sep="")]]))
  
  return(p1)
})

output$intPlateScatterPlot <- renderChart2({
  plot.data <- data()
  plot.data <- subset(plot.data, Plate==input$plateSelected & Experiment == input$experimentSelected)
  levels(plot.data$Control) <- c(levels(plot.data$Control), "Sample")
  plot.data[is.na(plot.data$Control), "Control"] <- "Sample"
  
  plot.data <- as.data.frame(plot.data)
  plot.data$Well.index <- seq(1, nrow(plot.data))
  plot.data <- plot.data[,c("Well.index", input$normalization, paste(input$normalization, "_sem", sep=""), "Control", "Accession", "Sample")]
  
  p <- highcharts.scatterplot.plate(plot.data, show.error=input$show.sem.in.hits)
  p$exporting(enabled=TRUE)
  return(p)
})

#create scatter plots of signal in all plates for each experiment using highcharts
#use taglist to bind them together
#TODO fix tooltip
observe({
  #if(!identical(input$mainNavbar, "Hits")) return(NULL)
  exp.data <- processedData()
  if(is.null(exp.data)) return(NULL)
  if(is.null(input$normalization)) sel.normalization <- "Raw"
  else sel.normalization <- input$normalization  
  
  for(experiment in unique(exp.data$Experiment))
  {
    local({
      plotData <- subset(exp.data, Experiment == experiment)
      plotName <- paste(experiment, "IntScatterPlot", sep="")
      output[[plotName]] <- renderChart({         
        formatted.plotData <-formatIntegerForPlot(plotData, "Plate")
        p <- hPlot(y= sel.normalization, x = "wellCount", data = formatted.plotData, type = "scatter", group = "Plate")
        #p$params$series[[1]]$data <- toJSONArray(formatted.plotData, json=F)
        #p$tooltip(formatter = "#! function(){ return this.point.x + ',k' + this.point.Sample + this.point.y;} !#")
        p$addParams(dom=plotName)      
        return(p)
      })
    })
  }
})

output$intScatterPlot <- renderChart({
  p <- hPlot(Raw ~ wellCount, data = processedData(), type = "scatter", group = "Plate")
  p$addParams(dom='intScatterPlot')
  return(p)
})

# Whole screen scatter plot, non interactive #
output$scatterPlot <- renderPlot({
  exp.data <- processedData()
  if(is.null(exp.data)) return(NULL)
  if(is.null(input$normalization)) sel.normalization <- "Raw"
  else sel.normalization <- input$normalization  
  q <- ggplot(exp.data, aes_string(x="wellCount", y=sel.normalization, color="Plate", shape="Replicate")) + geom_point() + geom_line(stat="hline", yintercept="mean", color="black", aes(group=interaction(Plate, Replicate))) 
  q <- q + theme_bw()
  q <- q + facet_grid(Readout~Experiment, scale="free")
  q <- q + guides(color=guide_legend(nrow=10, byrow=TRUE))
  print(q)
})

# Control plot for average plate values
output$plateMeanPlot <- renderPlot({
  exp.data <- processedData()
  if(is.null(exp.data)) return(NULL)
  plateMeanInfoFileName <- "help/plateMeanInfo.html"
  plateMeanInfoText <- readChar(plateMeanInfoFileName, file.info(plateMeanInfoFileName)$size)
  showshinyalert(session, "plateMeanInfo", plateMeanInfoText, "info")
  progress <- dummyProgressBar()
  on.exit(progress$close())
  
  q <- qplot(x=Plate, y=Raw, data=exp.data, geom="boxplot", color=Plate, shape=Replicate)
  q <- q + theme_bw() + facet_grid(Readout~Experiment, scale="free_y") + guides(color=guide_legend(nrow=10, byrow=TRUE))
  print(q)
})

# Plot for control performance (box plots) #
output$controlPlot <- renderPlot({  
  if(hasCtrls()) return(NULL)
  exp.data <- processedData()
  if(is.null(exp.data)) return(NULL)
  progress <- dummyProgressBar()
  on.exit(progress$close())
  plot.data <- filter(exp.data, !is.na(Control), tolower(Control) != "sample")
  q <- ggplot(plot.data, aes_string(x="Control", y="Raw", color="Plate", shape="Replicate")) + geom_boxplot()
  q <- q + theme_bw() + facet_grid(Readout~Experiment, scale="free_y") 
  q <- q + guides(color=guide_legend(nrow=10, byrow=TRUE))
  print(q)
})

# Plot for row and column effects over whole screen #
output$rowAndColumn <- renderPlot({
  exp.data <- processedData()
  if(is.null(exp.data)) stop("Please process the input data first")
  progress <- dummyProgressBar()
  on.exit(progress$close())
  p1 <- qplot(x=Row, y=centered, data=exp.data, geom="bar", stat="summary", fun.y="mean", fill=Replicate, position="dodge", main="Row Mean") + scale_y_continuous(labels=percent)
  p2 <- qplot(x=Column, y=centered, data=exp.data, geom="bar", stat="summary", fun.y="mean", fill=Replicate, position="dodge", main="Column Mean") + scale_y_continuous(labels=percent)
  p1 <- p1 + theme_bw() + facet_grid(Readout~Experiment) + scale_fill_brewer(palette="Accent")
  p2 <- p2 + theme_bw() + facet_grid(Readout~Experiment) + scale_fill_brewer(palette="Accent")
  print(grid.arrange(p1,p2))
})

KPM.modify.hits <- reactive({
  hits <- selectedHitList()
  
  if(input$accessionColType == "MI")
  {
    mimat <- as.data.frame(mirbaseMATURE)
    hits <- left_join(hits, mimat, by=c("mirna_id"))
  }
  return(hits)
})

#KPM miRNA target enrichment plot
#using igraph
output$KPM.plot.d3 <- renderForceNetwork({
  hits <- KPM.modify.hits()
  kpm.res <- KPM.result()
  if(is.null(kpm.res)) return(NULL)
  plot.miRNA.target.enrichment.graph.d3(kpm.res,targets.indicator.matrix(), hits, input$highlight.kpm_d3)  
})

#interactive using d3
output$KPM.plot.igraph <- renderPlot({
  hits <- KPM.modify.hits()
  kpm.res <- KPM.result()
  if(is.null(kpm.res)) return(NULL)
  plot.miRNA.target.enrichment.graph.igraph(kpm.res,targets.indicator.matrix(), hits)
})
