# find screening hits, e.g. the outliers #
outliers <- reactive({
  
  input$updateNormalization
  input$updateInclusion
  input$updateExclusion
  
  outl <- isolate({
    outl <- my.outliers(data(), input$method, input$margin)
    
    data <- data()
    outl[outl[[input$normalization]] > mean(data[[input$normalization]], na.rm=T),"category"] <- "promotor"
    outl[outl[[input$normalization]] < mean(data[[input$normalization]], na.rm=T),"category"] <- "suppressor"
    
    if(input$updateInclusion != 0 && nchar(input$include) > 0){
      if(length(grep(input$include, data$miRBase.ID.miRPlus.ID)) > 0){
        extra <- data[grep(input$include, data$miRBase.ID.miRPlus.ID),]
        extra$category <- "included"
        outl <- rbind(outl, extra)
      }
    }
    
    outl <- outl[,c(ncol(outl), seq(1:(ncol(outl)-1)))]
    
    return(outl)
  })
  return(outl)
})
