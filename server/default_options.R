#data option default
dataOptionDefaults <- reactive({
  if(input$dataset == "BCSC"){
    return(c("screenType" = "miRNA",
             "sampleCol" = "Sample", 
             "posColType" = "alpha",
             "posCol" = "Well.position",
             "accColType" = "MIMAT",
             "accCol" = "miRBase.accession",
             "measurementCol" = "CTB",
             "replicateCol" = "Replicate",
             "plateCol" = "Plate",
             "expCol" = "cellType",
             "ctrlCol" = "Control",
             "posCtrls" = "POS",
             "negCtrls" = "NEG",
             "rowCol" = "",
             "colCol" = "",
             "hasCtrls" = TRUE
    ))
  }
  
  else if(input$dataset == "A375_MTS"){
    return(c("screenType" = "miRNA",
             "sampleCol" = "miRNA", 
             "posColType" = "rowcol",
             "posCol" = "",
             "accColType" = "MI",
             "accCol" = "MI",
             "measurementCol" = "Raw",
             "replicateCol" = "Replicate",
             "plateCol" = "Plate",
             "expCol" = "Experiment",
             "ctrlCol" = "",
             "posCtrls" = "",
             "negCtrls" = "",
             "rowCol" = "Row",
             "colCol" = "Column",
             "hasCtrls" = FALSE
    ))
  }
  else return(NULL)
})
