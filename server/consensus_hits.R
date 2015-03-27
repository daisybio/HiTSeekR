# consensus hit list with all hits
outliers.all <- reactive({  
  exp.data <- data()
  #outl <- foreach(exp = unique(as.character(exp.data$Experiment)), .combine=rbind) %do%{
  exp.data <- subset(exp.data, Experiment==input$experimentSelected & Readout==input$readoutSelected)
  
  if(input$method=="Bayes" && is.null(negCtrl())) stop("negative control missing")
  
  progress <- shiny::Progress$new()
  progress$set(message = "Generating hit lists...", value = 0)
  on.exit(progress$close())
  
  #don't need to calculate default normalization twice
  normalizations <- input$multiNormalizations
  default.normalization <- input$normalization
  normalizations <- setdiff(normalizations, default.normalization)
  
  hits <- isolate({
    foreach(normalization = input$multiNormalizations, .combine=rbind) %do% 
    {      
      result <- find.hits.call(exp.data, input$method, input$margin, negCtrl(), normalization, NULL)
      result$method <- normalization
      return(result)
    }
  })
  
  #add default method hits
  if(default.normalization %in% input$multiNormalizations)
  {
    default.hits <- hit.detect()
    default.hits$method <- default.normalization
    hits <- rbind(hits, default.hits)
  }
  
  #count how often a sample is recognized as hit
  hits <- hits %>% count(Experiment, Readout, Plate, Well.position, Sample, method, category) 
  return(hits)
})

# reformat consensus hit list and apply threshold #
consensusHitList <- reactive({
  hits <- outliers.all()
    
  consensus <- hits %>% group_by(Experiment, Readout, Plate, Sample, Well.position) %>% 
    summarise(method=paste(method, collapse="/"), count=sum(n)) %>% 
    filter(count >= input$multiThreshold)
  
  consensus <- left_join(consensus, data(), by=c("Experiment", "Readout", "Plate", "Sample", "Well.position"))
  
  return(consensus)
})
