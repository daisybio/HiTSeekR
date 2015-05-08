ssmd <- function(exp.data, method, margin, neg.ctrl, signalColumn)
{    
  neg.ctrls <- filter(exp.data, Control %in% neg.ctrl)
  if(nrow(neg.ctrls) == 0) stop("At least one plate is missing negative controls")
  foreach(neg.ctrl = unique(neg.ctrls$Control), .combine=rbind) %do%
  {    
    neg.ctrl.data <- exp.data %>% filter(Control == neg.ctrl)
    if(nrow(neg.ctrl.data) < 3) stop("Cannot estimate variance in SSMD with < 3 negative control wells per plate")
    neg.ctrl.mean <- mean(neg.ctrl.data[[signalColumn]], na.rm=T)
    neg.ctrl.sd <- sd(neg.ctrl.data[[signalColumn]], na.rm=T)
    
    calc.ssmd <- function(exp.data)
    {
      sampleSignal <- exp.data[[signalColumn]]
      
      if(nrow(exp.data) == 1)
      {        
        return ((sampleSignal - neg.ctrl.mean) / (sqrt(2) * neg.ctrl.sd))
      }
      else
      {
        return ((mean(sampleSignal, na.rm=T) - neg.ctrl.mean) / sqrt((neg.ctrl.sd)^2 + sd(sampleSignal)^2))    
      }
    }
    
    exp.data <- exp.data %>% group_by(Sample, Well.position) %>% do(NEG.CTRL=neg.ctrl, SSMD=calc.ssmd(.))
    return(exp.data)
  }
}

find.hits.call <- function(exp.data, method, margin, neg.ctrl, signalColumn, updateProgress){
  
  if(method == "Bayes")
  {    
    #bayesian hit detection
    outl <- bayesianHitSelection(exp.data, neg.ctrl=neg.ctrl, signalColumn=signalColumn,alpha = 0.05, updateProgress=updateProgress)
  } 
  else if(method == "SSMD")
  {
    outl <- exp.data %>% group_by(Plate) %>% do(ssmd(., method, margin, neg.ctrl, signalColumn))
  }
  else{
    outl <- find.hits(exp.data, method, margin, signalColumn=signalColumn, updateProgress=updateProgress)    
    outl <- as.data.frame(outl)
  }  
  
  if(nrow(outl) == 0) stop("No hits were found with these settings.")
  outl[which(outl[,signalColumn] > mean(exp.data[,signalColumn], na.rm=T)),"category"] <- "promotor"
  outl[which(outl[,signalColumn] < mean(exp.data[,signalColumn], na.rm=T)),"category"] <- "suppressor"
  
  return(outl)
}

observeBayesAbuse <- observe({
  if(is.null(input$normalization) || is.null(input$method)) return(NULL)
  else if(input$normalization != "Raw" && input$method=="Bayes")
  {
    showshinyalert(session, "hits_error", "Bayes method should be used on raw data only", "danger")  
  } else if(input$method=="Bayes" && is.null(negCtrl())){
    showshinyalert(session, "hits_error", "Bayesian statistics are based on negative controls (used to calculate the priors). Select a negative control column in the DATA tab", "error")        
  }        
})

# find screening hits, e.g. the outliers #
hit.detect <- reactive({
  
  #input$updateNormalization
    
  if(input$method=="Bayes" && is.null(negCtrl())) stop("negative control missing")
  
  progress <- shiny::Progress$new()
  progress$set(message = "Finding hits...", value = 0)
  on.exit(progress$close())
  
  updateProgress <- function(value = NULL, detail = NULL) {
    if (is.null(value)) {
      value <- progress$getValue()
      value <- value + (progress$getMax() - value) / 5
    }
    progress$set(value = value, detail = detail)
  }  
    
  if(input$method != "SSMD") exp.data <- data()
  else exp.data <- processedData() #for SSMD we need the replicates
  #outl <- foreach(exp = unique(as.character(exp.data$Experiment)), .combine=rbind) %do%{
  exp.data <- subset(exp.data, Experiment==input$experimentSelected & Readout==input$readoutSelected)
  
  outl <- find.hits.call(exp.data, input$method, input$margin, negCtrl(), input$normalization, updateProgress)  
  
  return(outl)
})

outliers <- reactive({
  if(is.null(processedData())) return(NULL)
  outl <- hit.detect()
  
  #inclusion filter
  
  if(input$method == "Bayes"){
    exp.data <- outl
    outl <- outl[outl[[input$bayes_hypothesis]] > input$bayes_pval_cutoff,]  
  }
  else{
    exp.data <- data()
    exp.data <- subset(exp.data, Experiment==input$experimentSelected)
  }
  if(nchar(input$include) > 0){
    if(length(grep(input$include, exp.data$Sample)) > 0){
      extra <- exp.data[grep(input$include, exp.data$Sample),]
      extra$category <- "included"
      outl <- rbind(outl, extra)
    }
  }
  
  #exclusion filter
  
  if(nchar(input$exclude) > 0){
      outl <- outl[-grep(input$exclude, outl$Sample),]
  }

  
  #fix column order
  outl <- outl[,c(ncol(outl), seq(1:(ncol(outl)-1)))]
  
  return(outl)
})

#selected hit list in navbar_mirna_targets
selectedHitList <- reactive({
  if(input$useConsensus == "hit list") data <- outliers()
  else data <- consensusHitList()
  return(data)
})

htsanalyzerHitList <- reactive({
  if(input$htsanalyzer.useConsensus== "hit list") data <- outliers()
  else data <- consensusHitList()
  return(data)
})

#update select depending on what is selected in navbar_gsea
selectedHitListObserver <- observe({
input$htsanalyzer.useConsensus
  isolate({
    if(!is.null(input$useConsensus) && 
         input$htsanalyzer.useConsensus != input$useConsensus)
      updateSelectInput(session, "useConsensus", 
                      "Use normal hit list or consensus hit list for target identification?", 
                      c("hit list", "consensus hit list"), 
                      input$htsanalyzer.useConsensus)
  })
})

family.hitrate <- reactive({ 
 library(dplyr)
 all.data <- data()
 all.data <- all.data %>% dplyr::group_by(prefam_acc) %>% dplyr::summarise(library_count=dplyr::n_distinct(mirna_id)) 
 hits <- outliers()

 result <- formattedTable(hits, FALSE)
 result <- within(result, hits <- paste(category, Sample))
 result <- result %>% dplyr::group_by(prefam_acc, id) %>% dplyr::summarise(hits_count=n(),
                distinct_hits_count=n_distinct(mirna_id), 
                samples=paste(hits, collapse="")
 )
 
 final_result <- dplyr::left_join(result, all.data, by="prefam_acc")
 final_result <- final_result %>% mutate(family_coverage=distinct_hits_count/library_count)
 final_result <- final_result %>% filter(library_count > input$family_size_cutoff, family_coverage > (input$family_coverage_cutoff/100))
 if(nrow(final_result) == 0) stop("no families found")
 final_result <- as.data.frame(final_result)
 family_coverage <- as.numeric(final_result$family_coverage)
 final_result$family_coverage <- paste("<div style='background:#FDB462; text-align:center; border-radius: 15px; width:50px; height:25px;'>", 
                                       round(100*family_coverage, 0), "%",
                                       "</div>", sep="")
 final_result[family_coverage < 0.33, "family_coverage"] <-  sub("#FDB462", "#FF0000", final_result[family_coverage < 0.33, "family_coverage"])
 final_result[family_coverage > 0.66, "family_coverage"] <-  sub("#FDB462", "#40FF00", final_result[family_coverage > 0.66, "family_coverage"])
 final_result <- na.omit(final_result)
  
 return(as.data.frame(final_result))
})