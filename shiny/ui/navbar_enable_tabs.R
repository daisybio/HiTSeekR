observe({
  if(input$startButton > 0 && !is.null(processedData()))
  {    
    #enable quality control, hit list 
    session$sendCustomMessage(type = "enableNavTab", "2")
    session$sendCustomMessage(type = "enableNavTab", "3")
    session$sendCustomMessage(type = "enableNavTab", "4")
    
    #session$sendCustomMessage(type = "disableNavTab", "5")
    #session$sendCustomMessage(type = "disableNavTab", "6")
    #session$sendCustomMessage(type = "disableNavTab", "7")
    #session$sendCustomMessage(type = "disableNavTab", "8")
  }
})

observe({
  if(input$continueToQC)
  {
    updateTabsetPanel(session, "mainNavbar", selected = "Quality Control")
  }
})

observe({
  if(input$continueToNormalizationEffect){
    updateTabsetPanel(session, "mainNavbar", selected = "Normalization Effect")
  }
})

observe({
  if(input$continueToHitDiscovery){
    updateTabsetPanel(session, "mainNavbar", selected = "Hit Discovery")
  }
})

observe({
  if(input$continueToMiRNAs){
    updateTabsetPanel(session, "mainNavbar", selected = "microRNAs")
  }
})

observe({
  if(input$continueToDrugs){
    updateTabsetPanel(session, "mainNavbar", selected = "Small Compounds")
  }
})

observe({
  if(input$continueToGenes){
    updateTabsetPanel(session, "mainNavbar", selected = "Genes")
  }
})

observe({
  if(input$continueToGenes2){
    updateTabsetPanel(session, "mainNavbar", selected = "Genes")
  }
})

observe({
  if(input$continueToGenes3){
    updateTabsetPanel(session, "mainNavbar", selected = "Genes")
  }
})

observe({    
  if(input$mainNavbar == "Hit Discovery")
  { 
    #consensus hits based on hit discovery settings
    #session$sendCustomMessage(type = "enableNavTab", "5")
    
    #enable remaining tabs depending on screen type
    if(input$screenType == "miRNA"){
      session$sendCustomMessage(type = "enableNavTab", "5")  
    }
    else if(input$screenType == "compound"){
      session$sendCustomMessage(type = "enableNavTab", "6")  
    }   
    else{
      session$sendCustomMessage(type = "enableNavTab", "7")  
    }
  }
})

observe({
  if(input$mainNavbar == "microRNAs")
  {
    session$sendCustomMessage(type = "enableNavTab", "7")  
  }
  else return(NULL)
})

observe({
  if(input$mainNavbar == "Small Compounds")
  {
    session$sendCustomMessage(type = "enableNavTab", "7")  
  }
  else return(NULL)
})