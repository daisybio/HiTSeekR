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
    data <- as.data.frame(data)
    maxNumber <- max(as.integer(data[,columnToFormat]))
    data[,columnToFormat] <- sprintf(paste("%0", nchar(maxNumber), "d", sep=""),as.integer(data$Plate))
    data[,columnToFormat] <- as.factor(data[,columnToFormat])
    return(data)
}

# bar chart to visualize how many normalizations report a given hit #
# output$normcomparison <- renderChart({
#   exp.data <- outliers.all()
#   consensus <- consensusHitList()
#   
#   exp.data <- filter(exp.data, Sample %in% consensus$Sample)
#   p3 <- rPlot(x="Sample", y="n", color='method', data = exp.data, type = 'bar')
#   p3$addParams(dom='normcomparison', width="800", height="800")
#   return(p3)
# })

# Venn diagram to show overlap of hits in normalization methods #
# output$consensusVennDiagram <- renderPlot({
#   data <- outliers.all()
#   data <- split(data, data$method)
#   if(length(data) <= 5){
#     data <- lapply(data, function(x){na.omit(x$Sample)})
#     fill.colors <- c("purple", "cornflowerblue", "green", "yellow", "orange", "darkorchid1")
#     fill.colors <- fill.colors[1:length(data)]
#     p <- venn.diagram(data, filename = NULL, force.unique=T, col="transparent", 
#                       fill = fill.colors,
#                       alpha=0.5, fontfamily="serif")
#     grid.draw(p)
#   }
#   else textplt("this feature only supports up to 5 normalization methods")
# })

output$signalDistPlot <- renderPlot({
  plot.data <- processedData()
  progress <- dummyProgressBar()
  on.exit(progress$close())
  
  p <- ggplot(plot.data, aes_string(x=input$dataSelectedNormalization, fill="Replicate")) 
  p <- p + geom_density(alpha=.3)
  p <- p + theme_bw() + facet_grid(Experiment ~ Readout)
  return(p)  
}, height=600)

output$signalqqPlot <- renderPlot({
  plot.data <- processedData()
  progress <- dummyProgressBar()
  on.exit(progress$close())
  
  if(input$hasControls)
    opts <- aes_string(sample = input$dataSelectedNormalization, color="Control")
  else opts <- aes_string(sample = input$dataSelectedNormalization)
  p <- ggplot(plot.data, opts) + stat_qq()
  p <- p + theme_bw() + facet_grid(Experiment ~ Readout)
  
  return(p)  
}, height=600)

# Heatmap #
output$heatmapPlot <- renderPlot({  
  plot.data <- data()
  progress <- dummyProgressBar()
  on.exit(progress$close())
  
  ncol <- length(unique(plot.data$Plate))
  plot.data <- dplyr::filter(plot.data, Experiment == input$heatmapExperimentSelected, Readout == input$heatmapReadoutSelected)
  hits <- outliers()
  
  my.heatmap(plot.data, input$normalization, input$margin, input$method, hits, 
             showSampleLabels=input$showLabelsOnHeatmap, 
             signalColumn=input$normalization, 
             ncol=floor(sqrt(ncol)))
})

output$intHeatmapPlot <- renderChart2({
  plot.data <- processedData()
  plot.data <- subset(plot.data, Plate == input$plateSelected & 
                        Experiment == input$heatmapExperimentSelected  & 
                        Replicate == input$replicateSelected &
                        Readout == input$heatmapReadoutSelected)
  plot.data <- plot.data[,c("Column", "Row", input$dataSelectedNormalization, "Sample", "Accession")]
  colnames(plot.data)[1:3] <- c("x", "y", "value")
  plot.data <- arrange(plot.data, x, y)
  p <- highcharts.heatmap(as.data.frame(plot.data), input$dataSelectedNormalization, useWithShiny = T)  
  
  return(p)
})

renderChart3 <- function( expr, env = parent.frame(), quoted = FALSE ){
  func <- shiny::exprToFunction(expr, env, quoted)
  function() {
    rChart_ <- func()
    #cht_style <- sprintf("<style>.rChart {width: %spx; height: %spx} </style>", 
    #                     rChart_$params$width, rChart_$params$height)
    cht <- paste(
      capture.output(cat(
        rChart_$print()
        ,render_template(
          rChart_$templates$afterScript %||% 
            "<script></script>"
          , list(chartId = rChart_$params$dom, container = rChart_$container)
        )
        ,sep = ""
      ))
      , collapse = "\n")
    HTML(paste(cht, collapse = "\n"))
  }
}

# Interactive plots for screen hits #
output$scatterPlotHits <- renderChart3({
  outl <- outliers() %>% dplyr::select_("Sample", "Plate", "Well.position", "category", input$normalization, paste(input$normalization, "_sd", sep=""))
  outl <- unique(outl)

  p1 <- dPlot(y=input$normalization, x=c("Sample", "Plate", "Well.position"), z=paste(input$normalization, "_sd", sep=""), 
              data=outl, type="bubble", groups="category", height="600", width="1200",
              bounds = list(x=70, y=30, height="80%", width="90%"))  
  p1$addParams(dom='scatterPlotHits3')

  p1$legend(
    x = 580,
    y = 0,
    width = 50,
    height = 200,
    horizontalAlign = "left"
  )
  
  numOfSamples <- length(unique(as.character(outl$Sample)))
  if(numOfSamples > 500) stop("Too many hits to plot. Increase stringency.")
  if(numOfSamples > 20){
    numToSkip <- floor(numOfSamples / 20)
    axisFormatting <- paste(
    'var cleanAxis = function (axis, oneInEvery) {
            // This should have been called after draw, otherwise do nothing
            if (axis.shapes.length > 0) {
                // Leave the first label
                var del = 0;
                // If there is an interval set
                if (oneInEvery > 1) {
                    // Operate on all the axis text
                    axis.shapes.selectAll("text").each(function (d) {
                        // Remove all but the nth label
                        if (del % oneInEvery !== 0) {
                            this.remove();
                            // Find the corresponding tick line and remove
                            axis.shapes.selectAll("line").each(function (d2) {
                                if (d === d2) {
                                    this.remove();
                                }
                            });
                        }
                    del += 1;
                    });
                }
            }
        };
        
        cleanAxis(myChart.axes[0], ', numToSkip, ');', sep="")
  }
  else{
    axisFormatting <- ''
  } 
  
  p1$setTemplate(afterScript = 
                   paste('<script>
                          myChart.assignColor("included","#FDB462");
                          myChart.assignColor("promotor","#80B1D3");
                          myChart.assignColor("suppressor","#FB8072");',
                         axisFormatting,
    'myChart.draw();</script>', sep="")
  )

  return(p1)
})

output$intPlateScatterPlot <- renderChart2({
  plot.data <- data()
  progress <- dummyProgressBar()
  on.exit(progress$close())
  plot.data <- subset(plot.data, Plate==input$plateSelected & 
                        Experiment == input$heatmapExperimentSelected &
                        Readout == input$heatmapReadoutSelected)
  levels(plot.data$Control) <- c(levels(plot.data$Control), "Sample")
  plot.data[is.na(plot.data$Control), "Control"] <- "Sample"
  
  plot.data <- as.data.frame(plot.data)
  plot.data$Well.index <- seq(1, nrow(plot.data))
  plot.data <- plot.data[,c("Well.index", input$dataSelectedNormalization, paste(input$dataSelectedNormalization, "_sd", sep=""), "Control", "Accession", "Sample")]
  
  p <- highcharts.scatterplot.plate(plot.data, show.error=TRUE)
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
  progress <- dummyProgressBar()
  on.exit(progress$close())
  if(is.null(input$normalization)) sel.normalization <- "Raw"
  else sel.normalization <- input$normalization  
  
  for(experiment in unique(exp.data$Experiment))
  {
    for(readout in unique(exp.data$Readout))
    {
      local({
        plotData <- dplyr::filter(exp.data, Experiment == experiment, Readout == readout)
        plotName <- paste(experiment, readout, "IntScatterPlot", sep="")
        output[[plotName]] <- renderChart({         
          
          formatted.plotData <-formatIntegerForPlot(plotData, "Plate")
          p <- hPlot(y= sel.normalization, x = "sampleName", data = formatted.plotData, type = "scatter", group = "Plate")
          #p$params$series[[1]]$data <- toJSONArray(formatted.plotData, json=F)
          #p$tooltip(formatter = "#! function(){ return this.point.x + ',k' + this.point.Sample + this.point.y;} !#")
          p$addParams(dom=plotName)      
          return(p)
        })
      })
    }
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
  
  progress <- dummyProgressBar()
  on.exit(progress$close())
  
  if(is.null(input$dataSelectedNormalization)) sel.normalization <- "Raw"
  else sel.normalization <- input$dataSelectedNormalization
  exp.data$signal <- exp.data[[sel.normalization]]
  mean.data <- exp.data %>% 
    group_by(Experiment, Plate, Replicate) %>% 
    mutate(mean_signal = mean(signal, na.rm=T), minX = min(wellCount), maxX = max(wellCount))
  q <- ggplot(mean.data, aes_string(x="wellCount", y=sel.normalization, color="Plate", shape="Replicate")) + 
      geom_point() +
      geom_segment(aes(x = minX, xend = maxX, y = mean_signal, yend = mean_signal, group = interaction(Plate, Replicate)), colour = "black")
  
    #geom_line(stat="summary", fun.y="mean", color="black", aes(group=interaction(Plate, Replicate))) 
  q <- q + theme_bw()
  q <- q + facet_grid(Readout~Experiment, scales="free")
  q <- q + guides(color=guide_legend(nrow=10, byrow=TRUE))
  q <- q + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  q <- q + xlab("Well Count")    
  
  numOfPlates <- length(unique(exp.data$Plate))
  q <- q + scale_color_manual(values=rep(brewer.pal(8,"Dark2"), length.out=numOfPlates))
  print(q)
}, height=600)

getRainbowColors <- function(numOfCategories){
  #col <-rainbow(numOfCategories)
  col <- brewer.pal(8, "Dark2")
  col.index <- ifelse(seq(col) %% 2, 
                      seq(col), 
                      (seq(ceiling(length(col)/2), length.out=length(col)) %% length(col)) + 1)
  mixed <- col[col.index]
  return(mixed)
}

#control plot showing the separability of the positive and negative controls via SSMD
output$controlPerformancePlot <- renderPlot({
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  
  updateProgress <- function(value = NULL, detail = NULL) {
    if (is.null(value)) {
      value <- progress$getValue()
      value <- value + (progress$getMax() - value) / 5
    }
    progress$set(value = value, detail = detail)
  } 
  
  progress$set(message = "Generating plot...", value=0)
  
  exp.data <- processedData()    
  negCtrls <- negCtrl()
  posCtrls <- posCtrl()
  ctrls <- c(negCtrls, posCtrls)
  ctrl.data <- exp.data %>% filter(Control %in% ctrls)
  if(is.na(unique(ctrl.data$Sample))) ctrl.data$Sample <- ctrl.data$Control
  ctrl.data <- ungroup(ctrl.data)
  result <- ctrl.data %>% dplyr::group_by(Experiment, Readout, Plate, Replicate) 
  
  ssmdPlateCounter <<- 0
  ssmdPlateMax <<- length(dplyr::group_size(result))
  
  result <- result %>% do(ssmd(., negCtrls, "Raw", summarise.results=FALSE, updateProgress = updateProgress))    
  result <- as.data.frame(result)
  result <- result %>% filter(Sample != NEG.CTRL)
  
  p <- qplot(data=result, x=Plate, y=SSMD, size=I(3), color=Sample, shape=Replicate)  
  p <- p + facet_grid(Experiment + Readout ~ NEG.CTRL, scales="free") + theme_bw()
  p <- p + guides(color=guide_legend(nrow=10, byrow=TRUE))    
  p <- p + annotate("rect", xmin=0, xmax=length(unique(result$Plate))+1, ymin=-3, ymax=3, alpha=0.2, fill="red") 
  p <- p + annotate("rect", xmin=0, xmax=length(unique(result$Plate))+1, ymin=-3, ymax=-6, alpha=0.2, fill="orange") 
  p <- p + annotate("rect", xmin=0, xmax=length(unique(result$Plate))+1, ymin=3, ymax=6, alpha=0.2, fill="orange") 
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(p)
}, height=600)

#plot for correlation of replicates

correlationPlot <- function(exp.data, signalCol){
  exp.data$Replicate <- as.character(exp.data$Replicate)
  if(length(unique(exp.data$Replicate)) < 2) stop("No replicates have been found with the selected input settings.")
  
  exp.data$Signal <- exp.data[[signalCol]]
  combinations.of.replicates <- combn(unique(exp.data$Replicate),2)
  
  foreach(reps = iter(combinations.of.replicates, by="col")) %do%{
    
    repA <- dplyr::filter(exp.data, Replicate == reps[1]) %>% dplyr::select(Plate,Experiment, Readout, Well.position, Signal) %>% dplyr::rename(repA = Signal)
    repB <- dplyr::filter(exp.data, Replicate == reps[2]) %>% dplyr::select(Plate,Experiment, Readout, Well.position, Signal) %>% dplyr::rename(repB = Signal)
    
    plot.data <- dplyr::select(repA, Plate,Experiment, Readout, Well.position)
    plot.data <- dplyr::left_join(plot.data, repA, by=c("Plate", "Experiment", "Readout", "Well.position"))
    plot.data <- dplyr::left_join(plot.data, repB, by=c("Plate", "Experiment", "Readout", "Well.position"))
    
    p <- qplot(data=plot.data, x=repA, y=repB, xlab=reps[1], ylab=reps[2], color=Plate)
    p <- p + facet_grid(Experiment ~ Readout)
    p <- p + stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, aes(color=NULL), show_guide=F) 
    p <- p + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) 
    p <- p + geom_abline(intercept=0, slope=1, color="grey")    
    p <- p + guides(color=guide_legend(nrow=10, byrow=TRUE))+ theme_bw()
    
    numOfPlates <- length(unique(exp.data$Plate))
    p <- p + scale_color_manual(values=rep(getRainbowColors(8), length.out=numOfPlates))  
    return(p)
  }
}

#quality control tab
output$replicateCorrPlot <- renderPlot({
  exp.data <- processedData()
  progress <- dummyProgressBar()
  on.exit(progress$close())
  list.of.plots <- correlationPlot(exp.data, "Raw")
  do.call(grid.arrange, list.of.plots)
}, height=600)

#normalization tab
output$replicateCorrPlot2 <- renderPlot({
  exp.data <- processedData()
  progress <- dummyProgressBar()
  on.exit(progress$close())
  list.of.plots <- correlationPlot(exp.data, input$dataSelectedNormalization)
  do.call(grid.arrange, list.of.plots)
}, height=600)


# Control plot for plate signal spread (box plots)
plateSignalSpreadPlot <- function(exp.data, signal)
{
  q <- ggplot(exp.data) + geom_boxplot(aes_string(x="Plate", y=signal, color="Plate", shape="Replicate"))
  q <- q + theme_bw() + facet_grid(Readout~Experiment, scale="free_y") + guides(color=guide_legend(nrow=10, byrow=TRUE))
  q <- q + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  numOfPlates <- length(unique(exp.data$Plate))
  q <- q + scale_color_manual(values=rep(getRainbowColors(8), length.out=numOfPlates))  
  return(q)
}

#quality control tab
output$plateMeanPlot <- renderPlot({
  exp.data <- processedData()
  if(is.null(exp.data)) return(NULL)
  progress <- dummyProgressBar()
  on.exit(progress$close())
  
  q <- plateSignalSpreadPlot(exp.data, "Raw")
  print(q)
}, height=600)

#normalization tab
output$plateMeanPlot2 <- renderPlot({
  exp.data <- processedData()
  if(is.null(exp.data)) return(NULL)
  progress <- dummyProgressBar()
  on.exit(progress$close())
  
  q <- plateSignalSpreadPlot(exp.data, input$dataSelectedNormalization)
  print(q)
}, height=600)

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
  q <- q + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
  print(q)
}, height=600)

# Plot for row and column effects over whole screen #
output$rowAndColumn <- renderPlot({
  exp.data <- processedData()
  if(is.null(exp.data)) stop("Please process the input data first")
  progress <- dummyProgressBar()
  on.exit(progress$close())
  p1 <- ggplot(aes(x=Row, y=centered), data=exp.data) + 
    labs(title="Row Mean") + 
    geom_bar(stat="summary", aes(fill=Replicate), fun.y="mean", position="dodge") + 
    scale_y_continuous(labels=percent)
  p2 <- ggplot(aes(x=Column, y=centered), data=exp.data) +
    labs(title="Column Mean") +
    geom_bar(stat="summary", fun.y="mean", aes(fill=Replicate), position="dodge") + 
    scale_y_continuous(labels=percent)
  p1 <- p1 + theme_bw() + facet_grid(Readout~Experiment) + scale_fill_brewer(palette="Accent")
  p2 <- p2 + theme_bw() + facet_grid(Readout~Experiment) + scale_fill_brewer(palette="Accent")
  print(grid.arrange(p1,p2))
}, height=600)

#KPM miRNA target enrichment plot
#using d3
output$KPM.plot.d3 <- renderForceNetwork({
  graph.data <- kpm.graph.data()
  if(is.null(graph.data)) return(NULL)
  plot.kpm.d3(graph.data, KPM.modify.hits(), screenType=input$screenType)
})

#interactive using igraph
output$KPM.plot.igraph <- renderPlot({
  graph.data <- kpm.graph.data()
  if(is.null(graph.data)) return(NULL)
  plot.kpm.igraph(graph.data, KPM.modify.hits(), input$screenType)
})
