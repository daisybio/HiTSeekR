# find screening hits, e.g. the outliers #
outliers <- reactive({
  
  input$updateNormalization
  input$updateInclusion
  input$updateExclusion
  
  outl <- isolate({
    exp.data <- data()
    #outl <- foreach(exp = unique(as.character(exp.data$Experiment)), .combine=rbind) %do%{
    outl <- my.outliers(subset(exp.data, Experiment==input$experimentSelected), input$method, input$margin, signalColumn=input$normalization)    
         
    if(nrow(outl) == 0) return(outl)
    outl[which(outl[,input$normalization] > mean(exp.data[,input$normalization], na.rm=T)),"category"] <- "promotor"
    outl[which(outl[,input$normalization] < mean(exp.data[,input$normalization], na.rm=T)),"category"] <- "suppressor"
    
    if(input$updateInclusion != 0 && nchar(input$include) > 0){
      if(length(grep(input$include, data$Sample)) > 0){
        extra <- data[grep(input$include, data$Sample),]
        extra$category <- "included"
        outl <- rbind(outl, extra)
      }
    }
    
    outl <- outl[,c(ncol(outl), seq(1:(ncol(outl)-1)))]
    
    return(outl)
  })
  return(outl)
})

family.hitrate <- reactive({
 library(mirbase.db)
 all.data <- data()
 all.data <- all.data %>% group_by(prefam_acc) %>% summarise(library_count=n_distinct(mirna_id)) 
 hits <- outliers()

 result <- formattedTable()
 result <- within(result, hits <- paste(category, mature_name))
 result <- result %>% group_by(prefam_acc, id) %>% summarise(hits_count=n(),
                distinct_hits_count=n_distinct(mirna_id), 
                category_mature_mirna=paste(hits, collapse="")
 )
 #
 final_result <- left_join(result, all.data, by="prefam_acc")
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
 
 showshinyalert(session, "miRNA_family_info", "library_count refers to ", "info")
 
 return(as.data.frame(final_result))
})