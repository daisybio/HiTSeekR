# find screening hits, e.g. the outliers #
outliers <- reactive({
  
  input$updateNormalization
  input$updateInclusion
  input$updateExclusion
  
  outl <- isolate({
    exp.data <- data()
    
    outl <- my.outliers(exp.data, input$method, input$margin, signalColumn=input$normalization)    
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
