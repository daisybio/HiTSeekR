#data option default
dataOptionDefaults <- reactive({  
  if(datasetName() == "BCSC"){
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
             "hasCtrls" = "TRUE"
    ))
  }  
  else if(datasetName() == "A375_MTS"){
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
             "hasCtrls" = "FALSE"
    ))
  }
  else if(datasetName() == "DM_Kc167"){
    return(c("screenType" = "siRNA",
             "sampleCol" = "HFAID", 
             "posColType" = "alpha",
             "posCol" = "well",
             "accColType" = "FlybaseCG",
             "accCol" = "GeneID",
             "measurementCol" = "value",
             "replicateCol" = "replicate",
             "plateCol" = "plate",
             "expCol" = "experiment",
             "ctrlCol" = "controlStatus",
             "posCtrls" = "pos",
             "negCtrls" = "neg",
             "rowCol" = "",
             "colCol" = "",
             "hasCtrls" = "TRUE"
    ))
  }
  else return(c("screenType" = "",
                    "sampleCol" = "", 
                    "posColType" = "",
                    "posCol" = "",
                    "accColType" = "",
                    "accCol" = "",
                    "measurementCol" = "",
                    "replicateCol" = "",
                    "plateCol" = "",
                    "expCol" = "",
                    "ctrlCol" = "",
                    "posCtrls" = "",
                    "negCtrls" = "",
                    "rowCol" = "",
                    "colCol" = "",
                    "hasCtrls" = "FALSE"
  ))
})
