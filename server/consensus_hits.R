# consensus hit list with all hits
outliers.all <- reactive({
  hits <- data.frame()
  my.data <- data()
  input$updateNormalization
  
  for(method in input$multiNormalizations)
  {
    my.data$signal <- my.data[[method]]
    outl <- isolate(my.outliers(my.data, input$method, input$margin, withControls=input$includeControls))
    outl$method <- method
    
    outl[outl[[method]] > mean(my.data[[method]], na.rm=T),"category"] <- "promotor"
    outl[outl[[method]] < mean(my.data[[method]], na.rm=T),"category"] <- "suppressor"
    
    hits <- rbind(hits, outl)
  }
  
  return(ddply(hits, .(Sample, method, category), summarise, count=1))
})

# reformat consensus hit list and apply threshold #
consensusHitList <- reactive({
  consensus <- ddply(outliers.all(), .(Sample, category), summarise, method=paste(method, collapse="/"), count=length(count))
  consensus <- subset(consensus, count >= input$multiThreshold)
  consensus <- merge(consensus, data(), by="Sample")
  consensus <- subset(consensus, select =c(2, 1, seq(3,(ncol(consensus)))))
  consensus$signal <- NULL
  consensus$signal.sem <- NULL
  return(consensus)
})
