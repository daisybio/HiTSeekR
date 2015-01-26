#format integers for plots
formatIntegerForPlot <- function(data, columnToFormat){    
    maxNumber <- max(as.integer(data[,columnToFormat]))
    data[,columnToFormat] <- sprintf(paste("%0", nchar(maxNumber), "d", sep=""),as.integer(data$Plate))
    data[,columnToFormat] <- as.factor(data[,columnToFormat])
    return(data)
}

# bar chart to visualize how many normalizations report a given hit #
output$normcomparison <- renderChart({
  data <- outliers.all()
  data <- subset(data, Sample %in% consensusHitList()$Sample)
  p3 <- rPlot(x="Sample", y="count", color='method', data = data, type = 'bar')
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
  p1 <- dPlot(y=input$normalization, x=c("Sample", "Well.position"), z=paste(input$normalization, ".sem", sep=""), 
              data=outl, type="bubble", groups="category", height="800", width="100%",
              bounds = list(x=70, y=30, height="600", width="90%"))  
  p1$addParams(dom='scatterPlotHits')
  if(input$show.sem.in.hits) p1$zAxis(type="addMeasureAxis", overrideMax = 2*max(outl[[paste(input$normalization, ".sem", sep="")]]))
  
  return(p1)
})

output$intPlateScatterPlot <- renderChart2({
  plot.data <- data()
  plot.data <- subset(plot.data, Plate==input$plateSelected & Experiment == input$experimentSelected)
  levels(plot.data$Control) <- c(levels(plot.data$Control), "Sample")
  plot.data[is.na(plot.data$Control), "Control"] <- "Sample"
  
  plot.data <- as.data.frame(plot.data)
  plot.data$Well.index <- seq(1, nrow(plot.data))
  plot.data <- plot.data[,c("Well.index", input$normalization, paste(input$normalization, ".sem", sep=""), "Control", "Accession", "Sample")]
  
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
### DEPRECATED ###
output$scatterPlot <- renderPlot({
  data <- processedData()
  q <- ggplot(data, aes_string(x="wellCount", y=input$normalization, color="Plate", shape="Replicate")) + geom_point() + geom_line(stat="hline", yintercept="mean", color="black", aes(group=interaction(Plate, Replicate))) 
  q <- q + theme_bw()
  if("Experiment" %in% colnames(data)) q <- q + facet_wrap(~Experiment, ncol=1, scales="free")
  print(q)
})

# Plot for control performance (box plots) #
output$controlPlot <- renderPlot({
  data <- processedData()
  q <- ggplot(subset(data, !is.na(Control)), aes_string(x="Sample", y=input$normalization, color="Plate", shape="Replicate")) + geom_boxplot()
  
  print(q)
})

# Plot for row and column effects over whole screen #
output$rowAndColumn <- renderPlot({
  p1 <- qplot(x=Row, y=centered, data=processedData(), geom="bar", stat="summary", fun.y="mean", fill=Replicate, position="dodge", main="Row Mean") + scale_y_continuous(labels=percent)
  p2 <- qplot(x=Column, y=centered, data=processedData(), geom="bar", stat="summary", fun.y="mean", fill=Replicate, position="dodge", main="Column Mean") + scale_y_continuous(labels=percent)
  print(grid.arrange(p1,p2))
})

# Plot for comparison of normalization methods in plate viewer #
output$normalizationComparison <- renderPlot({
  data <- processedData()
  data$signal <- NULL
  p2 <- qplot(x=wellCount, y=value, data=melt(subset(data, Plate==input$plateSelected), id.vars=c("Plate", "Well.position", "wellCount", "Replicate", "Sample", "Accession", "Sample", "Control", "Row", "Column")), color=Replicate, main="Comparison of different normalization methods") + facet_wrap(~variable, scales="free", ncol=3) + geom_smooth(method="loess")
  print(p2)
})

#KPM miRNA target enrichment plot
output$KPM.plot <- renderPlot({

  hits <- selectedHitList()
  
  if(input$accessionColType == "MI")
  {
    mimat <- as.data.frame(mirbaseMATURE)
    hits <- left_join(hits, mimat, by=c("mirna_id"))
  }
  plot.miRNA.target.enrichment.graph(KPM.result(),targets.indicator.matrix(), hits)
})
