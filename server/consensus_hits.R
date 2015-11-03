# consensus hit list with all hits
outliers.all <- reactive({  
  exp.data <- data()  
  rep.data <- processedData()
  
  if(input$method=="Bayes" && is.null(negCtrl())) stop("negative control missing")
  
  progress <- shiny::Progress$new()
  progress$set(message = "Generating hit lists...", value = 0)
  on.exit(progress$close())
  
  #don't need to calculate default normalization twice
  normalizations <- input$multiNormalizations
  default.normalization <- input$normalization
  normalizations <- setdiff(normalizations, default.normalization)
  
  hits <- isolate({
    foreach(exp = input$experimentSelected, .combine=rbind) %do%{
      foreach(rdt = input$readoutSelected, .combine=rbind) %do% {
        foreach(normalization = input$multiNormalizations, .combine=rbind) %do% 
        {      
          exp.d <- dplyr::filter(exp.data, Experiment == exp, Readout == rdt)
          rep.d <- dplyr::filter(rep.data, Experiment == exp, Readout == rdt)
          result <- find.hits.call(exp.d, rep.d, input$method, input$margin, negCtrl(), normalization, NULL)
          result$method <- normalization
          return(result)
        }
      }
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
  hits <- hits %>% dplyr::count(Experiment, Readout, Plate, Well.position, Sample, method, category) 
  return(hits)
})

# reformat consensus hit list and apply threshold #
consensusHitList <- reactive({
  hits <- outliers.all()
    
  consensus <- hits %>% dplyr::group_by(Experiment, Readout, Plate, Sample, Well.position) %>% 
    dplyr::summarise(method=paste(method, collapse="/"), count=sum(n)) %>% 
    filter(count >= input$multiThreshold)
  exp.data <- data()
  consensus <- dplyr::left_join(consensus, exp.data, by=c("Experiment", "Readout", "Plate", "Sample", "Well.position"))
  
  return(consensus)
})
