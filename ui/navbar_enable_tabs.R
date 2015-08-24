observe({
  if(input$startButton > 0)
  {    
    #enable quality control, hit list 
    session$sendCustomMessage(type = "enableNavTab", "2")
    session$sendCustomMessage(type = "enableNavTab", "3")
    session$sendCustomMessage(type = "enableNavTab", "4")
  }
})

observe({    
  if(input$mainNavbar == "Hit Discovery")
  { 
    #consensus hits based on hit discovery settings
    session$sendCustomMessage(type = "enableNavTab", "5")
    
    #enable remaining tabs depending on screen type
    if(input$screenType == "miRNA"){
      session$sendCustomMessage(type = "enableNavTab", "6")  
    }
    else if(input$screenType == "compound"){
      session$sendCustomMessage(type = "enableNavTab", "7")  
    }   
    else{
      session$sendCustomMessage(type = "enableNavTab", "8")  
    }
  }
})

observe({
  if(input$mainNavbar == "microRNAs")
  {
    session$sendCustomMessage(type = "enableNavTab", "8")  
  }
  else return(NULL)
})

observe({
  if(input$mainNavbar == "Small Compounds")
  {
    session$sendCustomMessage(type = "enableNavTab", "8")  
  }
  else return(NULL)
})