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
  input$updateNormalization
  isolate(my.heatmap(data(), input$normalization, input$margin, input$method, outliers(), colorA=input$colorA, colorB=input$colorB, showSampleLabels=input$showLabelsOnHeatmap))
})

# Interactive plots for screen hits #
output$scatterPlotHits <- renderChart({
  outl <- outliers()
  p1 <- dPlot(y=input$normalization, x=c("Sample", "Well.position"), z=paste(input$normalization, ".sem", sep=""), data=outl, type="bubble", groups="category", width=1000, height=600)  
  p1$addParams(dom='scatterPlotHits')
  if(input$showSEM) p1$zAxis(type="addMeasureAxis", overrideMax = 2*max(outl[[paste(input$normalization, ".sem", sep="")]]))
  
  return(p1)
})

# Plate viewer interactive scatter plot #
output$scatterPlotI <- renderChart({
  data <- subset(data(), Plate==input$plateSelected)
  p2 <- dPlot(y=input$normalization, x=c("Sample", "Well.position"), z=paste(input$normalization, ".sem", sep=""), data=data, type="bubble", groups="Control", width=1000, height=600)
  p2$xAxis(type="addCategoryAxis", orderRule="wellCount")
  #p2$yAxis(type="addMeasureAxis")
  p2$addParams(dom='scatterPlotI')
  if(input$showSEM) p2$zAxis(type="addMeasureAxis", overrideMax = 2*max(data[[paste(input$normalization, ".sem", sep="")]]))
  return(p2)
})

# Whole screen scatter plot #
output$scatterPlot <- renderPlot({
  data <- processedData()
  q <- ggplot(data, aes_string(x="wellCount", y=input$normalization, color="Plate", shape="Replicate")) + geom_point() + geom_line(stat="hline", yintercept="mean", color="black", aes(group=interaction(Plate, Replicate)))
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
  p2 <- qplot(x=wellCount, y=value, data=melt(subset(data, Plate==input$plateSelected), id.vars=c("Plate", "Well.position", "wellCount", "Replicate", "miRBase.ID.miRPlus.ID", "miRBase.accession", "Sample", "Control", "Row", "Column")), color=Replicate, main="Comparison of different normalization methods") + facet_wrap(~variable, scales="free", ncol=3) + geom_smooth(method="loess")
  print(p2)
})
