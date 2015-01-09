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
             "negCtrls" = "NEG"
    ))
  }
  else return(NULL)
})
