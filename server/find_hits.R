#reactive element that post-processes the hit list
outliers <- reactive({
  
  #first check if a hit list was uploaded directly
  if(input$isHitList){
    outl <- data()
    outl$category <- "included"
    return(outl)
  }
  
  #nothing to process, return null
  if(is.null(processedData())) return(NULL)
  
  #here, the actual hit detection is done
  outl <- hit.detect()
  
  #nothing to show so return null
  if(is.null(outl)) return(NULL)
  
  #filtering
  if(input$method == "Bayes"){
    if(input$normalization != "Raw") return(NULL)
    #depending on the effect we are interested in we have to use a different column for filtering
    #each effect is represented as its own bayesian hypothesis in a column
    bayes_hypothesis <- switch(input$effect,
                               "effect" = "p_effect",
                               "suppressor" =  "p_suppressor",
                               "promotor" = "p_promotor")
    
    #keep the p-values of all results here for the subsequent inclusion filter (regular expression)
    exp.data <- outl
    #extract the hits using the selected p-value cutoffs
    outl <- outl[which(outl[[bayes_hypothesis]] > input$bayes_pval_cutoff),]  
  }
  else{ #filter all other methods for effect category and prepare exp.data based on all sample values
    exp.data <- data()    
    exp.data <- dplyr::filter(exp.data, Experiment %in% input$experimentSelected, 
                              Readout %in% input$readoutSelected)
    if(input$effect != "effect") outl <- dplyr::filter(outl, category == input$effect)
  }
  
  #inclusion filter based on regular expression. samples are added to hit list with category included
  if(nchar(input$include) > 0){
    if(length(grep(input$include, exp.data$Sample)) > 0){
      extra <- exp.data[grep(input$include, exp.data$Sample),]
      extra$category <- "included"
      outl <- rbind(outl, extra)
    }
  }
  
  #exclusion filter using regular expression
  if(nchar(input$exclude) > 0){
    outl <- outl[-grep(input$exclude, outl$Sample),]
  }
  
  #fix column order
  outl <- outl[,c(ncol(outl), seq(1:(ncol(outl)-1)))]
  
  #sort after signal
  outl <- arrange_(outl, input$normalization)  
  
  #differential screening. remove hits found in the reference screen
  if(input$differentialScreening){
    if(length(input$experimentSelected) == 1 && length(input$readoutSelected) == 1){
      showshinyalert(session, "hits_error", "You have to select at least two different readouts or experiments to apply differential screening.", "danger")
      return(NULL)
    }
    else{      
      refSet <- dplyr::filter(outl, Experiment == input$referenceExperiment, Readout == input$referenceReadout)
      outl <- dplyr::anti_join(outl, refSet, by=c("Plate", "Well.position"))
    }
  }
  
  #return results
  return(as.data.frame(outl))
})

#Observer to check if Bayes is used correctly, e.g. on raw data.
observeBayesAbuse <- observe({
  if(is.null(input$normalization) || is.null(input$method)) return(NULL)
  else if(input$normalization != "Raw" && input$method=="Bayes")
  {
    updateSelectInput(session, "normalization", "Normalization", normalizationChoices(), "Raw")
    showshinyalert(session, "hits_error", "Bayes method should be used on raw data only", "danger")  
  } else if(input$method=="Bayes" && is.null(negCtrl())){
    showshinyalert(session, "hits_error", "Bayesian statistics are based on negative controls (used to calculate the priors). Select a negative control column in the DATA tab", "danger")        
  }        
})

#find screening hits using the selected method and parameters
hit.detect <- reactive({
  #We don't want Bayes method to start just by selecting it because it's computationally expensive.
  #We thus check if a dedicated apply button has been pressed since this code was executed last.
  #And we show a warning to make the user aware that he or she has to press a button.
  if(!exists("bayesButtonCounter")) bayesButtonCounter <<- 0
  if(input$method == "Bayes" && input$computeBayes == bayesButtonCounter){
    showshinyalert(session, "hits_error", "Computation of Bayesian statistics is computationally expensive and will take a while to compute. Press 'Apply Bayes method' to trigger the computation.", "danger")       
    return(NULL)
  } 
  else{
    bayesButtonCounter <<- input$computeBayes
  }
  
  #create a progress bar
  progress <- shiny::Progress$new()
  progress$set(message = "Hit discovery in progress...", value = 0)
  on.exit(progress$close())
  
  updateProgress <- function(value = NULL, detail = NULL) {
    if (is.null(value)) {
      value <- progress$getValue()
      value <- value + (progress$getMax() - value) / 5
    }
    progress$set(value = value, detail = detail)
  } 

  #gather input data
  merged.data <- data()
  exp.data <- processedData() 
  
  #to trigger reactive update
  margin <- input$margin
  
  #split input data by experiment and readout and process each individually
  
  errors.log <<- ""
  
  outl <- foreach(exp = input$experimentSelected, .combine=rbind) %do%{
    foreach(rdt = input$readoutSelected, .combine=rbind) %do% {
      
      #filter for specific screen
      it.data <- dplyr::filter(exp.data, Readout==rdt, Experiment==exp)
      m.data <- dplyr::filter(merged.data, Readout==rdt, Experiment==exp)   
      
      if(!is.null(input$referenceExperiment) && !is.null(input$referenceReadout) && input$differentialScreening){            
      
        #in case of a differential screening setup we might want to use a different detection margin for 
        #target and reference screen
        if(exp == input$referenceExperiment && rdt == input$referenceReadout){
          margin <- input$diffMargin   
        }
        else margin <- input$margin
      }
      
      #function to apply hit detection
      result <- NULL
      
      tryCatch({
        result <- find.hits.call(m.data, it.data, input$method, margin, negCtrl(), input$normalization, updateProgress, upperCutoff=input$upperCutoff, lowerCutoff=input$lowerCutoff)
      },
      error = function(e) { 
          errors.log <<- paste(errors.log, "For experiment ", exp, " and readout ", rdt, e$message, "<br/>")        
        }
      )
    }
      return(result)
  }
  
  #check how many of the hits are not mapped to an unambigious identifier
  if(input$screenType == 'siRNA')
  {
    na.count <- length(which(is.na(outl$gene_id)))
    possible_reason <- "Gene symbols not accepted by HUGO are sometimes ambigious or do not match to an entrez gene id. Typical examples are hypothetical genes starting with LOC."
  }
  else if(input$screenType == 'miRNA')
  {
    na.count <- length(which(is.na(outl$mirna_id)))
    possible_reason <- "miRNA names other than MI or MIMAT are often outdated. This is not necessarily a problem, but for some miR identifiers it is not clear whether the 3p or 5p strand is meant. Moreover, some hits are ignored because they belong to dead miRNA entries."
  }
  else if(input$screenType == 'compound')
  {
    #convert ids to pubchem compound ids (CID) as used in STITCH
    hits <- convertToCid(outl, input$accessionColType)
    
    #add CID prefix and leading zeros for STITCH db
    hits[!is.na(hits$PubChem_CID), "PubChem_CID"] <- paste("CID", formatC(as.integer(hits[!is.na(hits$PubChem_CID), "PubChem_CID"]), width=9, flag="0"), sep="")
    
    outl <- hits[,c(1:(ncol(hits)-2), ncol(hits), ncol(hits) -1)]
    na.count <- length(which(is.na(outl$PubChem_CID)))
    
    possible_reason <- "Some of the identifiers used here have no known match to a PubChem compound id (CID)."
  }
  
  if(na.count > 0){
    errors.log <-  paste(errors.log, "Warning: ", na.count, " of the selected hits could not be mapped to a suitable identifier and will thus be ignored in subsequent analyses. Possible reason: ", possible_reason, sep="")
  }
  
  #output errors and warnings
  if(errors.log != "") showshinyalert(session, "hits_error", errors.log, "danger")
  else hideshinyalert(session, "hits_error")
  
  #return results    
  return(outl)
})