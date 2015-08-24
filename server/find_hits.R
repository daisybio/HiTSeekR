ssmd <- function(exp.data, neg.ctrl, signalColumn, summarise.results=TRUE)
{        
  neg.ctrls <- dplyr::filter(exp.data, Control %in% neg.ctrl)
  if(nrow(neg.ctrls) == 0) stop("At least one plate is missing negative controls")
  result <- foreach(nc = neg.ctrl, .combine=rbind) %do%
  {      
    neg.ctrl.data <- exp.data %>% dplyr::filter(Control == nc)
    if(nrow(neg.ctrl.data) < 3) stop("Cannot estimate variance in SSMD with < 3 negative control wells per plate")
    neg.ctrl.mean <- mean(neg.ctrl.data[[signalColumn]], na.rm=T)
    neg.ctrl.sd <- sd(neg.ctrl.data[[signalColumn]], na.rm=T)
    
    calc.ssmd <- function(x)
    {      
      sampleSignal <- x[[signalColumn]]
      
      if(nrow(x) == 1)
      {        
        return ((sampleSignal - neg.ctrl.mean) / (sqrt(2) * neg.ctrl.sd))
      }
      else
      {
        return ((mean(sampleSignal, na.rm=T) - neg.ctrl.mean) / sqrt((neg.ctrl.sd)^2 + sd(sampleSignal, na.rm=T)^2))    
      }
    }    
    result <- as.data.frame(exp.data)
    result <- result %>% dplyr::group_by(Sample) %>% 
      dplyr::do(NEG.CTRL=nc, SSMD=calc.ssmd(.))       
    return(na.omit(result))
  }
  if(summarise.results)
    result <- result %>% dplyr::group_by(Sample) %>% dplyr::summarise(SSMD = mean(unlist(SSMD), na.rm=T))
  else{
    result <- result %>% dplyr::mutate(SSMD = unlist(SSMD), NEG.CTRL = unlist(NEG.CTRL))
  }
  return(as.data.frame(result))
}

find.hits.call <- function(exp.data, rep.data, method, margin, neg.ctrl, signalColumn, updateProgress){
  
  if(method == "Bayes")
  {    
    #bayesian hit detection
    outl <- bayesianHitSelection(exp.data, neg.ctrl=neg.ctrl, signalColumn=signalColumn,alpha = 0.05, updateProgress=updateProgress)
  } 
  else if(margin == 0) outl <- exp.data
  else if(method == "SSMD")
  {        
    result <- rep.data %>% dplyr::group_by(Plate) %>% do(ssmd(., neg.ctrl, signalColumn))
    outl <- exp.data
    outl <- dplyr::left_join(outl,result, by=c("Plate", "Sample"))
    outl <- outl %>% dplyr::filter(abs(SSMD) >= margin)    
    outl <- as.data.frame(outl)
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

  merged.data <- data()
  exp.data <- processedData() 
  outl <- foreach(exp = input$experimentSelected, .combine=rbind) %do%{
    foreach(rdt = input$readoutSelected, .combine=rbind) %do% {
      it.data <- dplyr::filter(exp.data, Readout==rdt, Experiment==exp)
      m.data <- dplyr::filter(merged.data, Readout==rdt, Experiment==exp)
      
      if(!is.null(input$referenceExperiment) && !is.null(input$referenceReadout)){            
        if(exp == input$referenceExperiment && rdt == input$referenceReadout){
          margin <- input$diffMargin   
        }
        else margin <- input$margin
      }      
      find.hits.call(m.data, it.data, input$method, margin, negCtrl(), input$normalization, updateProgress)      
    }
  }  
    
  return(outl)
})

outliers <- reactive({
  if(is.null(processedData())) return(NULL)
  outl <- hit.detect()
  
  #inclusion filter
  
  if(input$method == "Bayes"){
    bayes_hypothesis <- switch(input$effect,
           "effect" = "p_effect",
           "suppressor" =  "p_suppressor",
           "promotor" = "p_promotor")
    exp.data <- outl
    outl <- outl[outl[[bayes_hypothesis]] > input$bayes_pval_cutoff,]  
  }
  else{
    exp.data <- data()    
    exp.data <- dplyr::filter(exp.data, Experiment %in% input$experimentSelected, 
                              Readout %in% input$readoutSelected)
    if(input$effect != "effect") outl <- dplyr::filter(outl, category == input$effect)
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
  
  #sort after signal
  outl <- arrange_(outl, input$normalization)  
  
  #differential screening
  if(input$differentialScreening){
    if(length(input$experimentSelected) == 1 && length(input$readoutSelected) == 1)
      stop("You have to select at least two different readouts or experiments to apply differential screening.")
    else{      
      refSet <- dplyr::filter(outl, Experiment == input$referenceExperiment, Readout == input$referenceReadout)
      outl <- dplyr::anti_join(outl, refSet, by=c("Plate", "Well.position"))
    }
  }
  
  return(as.data.frame(outl))
})

genes.indicator.matrix <- reactive({
  if(input$KPM.useConsensus == "hit list") hits <- outliers()
  else hits <- consensusHitList()
  all.samples <- data()
  
  gene_ids <- na.omit(unique(all.samples$gene_id))
  ind.matrix <- as.matrix(rep(0, length(gene_ids)))
  row.names(ind.matrix) <- gene_ids
  ind.matrix[which(gene_ids %in% hits$gene_id), 1] <- 1
  
  return(ind.matrix)
})

htsanalyzerHitList <- reactive({
  if(input$htsanalyzer.useConsensus== "hit list") data <- outliers()
  else data <- consensusHitList()
  return(data)
})

KPM.selectedHitList <- reactive({
  if(input$KPM.useConsensus== "hit list") data <- outliers()
  else data <- consensusHitList()
  return(data)  
})

#update select depending on what is selected in navbar_gsea
# selectedHitListObserver <- observe({
# input$htsanalyzer.useConsensus
#   isolate({
#     if(!is.null(input$useConsensus) && 
#          input$htsanalyzer.useConsensus != input$useConsensus)
#       updateSelectInput(session, "useConsensus", 
#                       "Use normal hit list or consensus hit list for target identification?", 
#                       c("hit list", "consensus hit list"), 
#                       input$htsanalyzer.useConsensus)
#   })
# })

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