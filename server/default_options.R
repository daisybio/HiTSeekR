#data option default
dataOptionDefaults <- reactive({     
  
  if(datasetName() == "TNFa_Casp4"){
    return(list("screenType" = "siRNA",
                "sampleCol" = "Target.gene.symbol", 
                "posColType" = "numeric",
                "posCol" = "Position",
                "accColType" = "GeneSymbol",
                "accCol" = "Target.gene.symbol",
                "measurementCol" = "value",  
                "replicateCol" = "replicate",
                "plateCol" = "Plate.ID",
                "expCol" = "channel",
                "ctrlCol" = "Silencing.RNA.reagent.identifier",
                "posCtrls" = c("rela", "tnfr1"),
                "negCtrls" = c("sicon", "lrp5"),
                "rowCol" = "",
                "colCol" = "",
                "hasCtrls" = "TRUE"
    ))
  }
  else if(datasetName() == "HCC_vorinostat_miRNA"){
    return(list("screenType" = "miRNA",
             "sampleCol" = "X2.miRNAs.information.STRING....", 
             "posColType" = "alpha",
             "posCol" = "X15.WELL_LOCATION.STRING....",
             "accColType" = "MIMAT",
             "accCol" = "MIMAT",
             "measurementCol" = c(#"X5.AVERAGE_CTF_ROBUST_Z_SCORE.FLOAT....", 
                                  #"X6.AVERAGE_CTF_NORMALISED_TO_MOCK.FLOAT....", 
                                  "X7.AVERAGE_CTF_FLUORESCENCE_UNITS.INT....",
                                  #"X8.AVERAGE_CASP_ROBUST_Z_SCORE.FLOAT....",
                                  #"X9.AVERAGE_CASP_NORMALISED_TO_MOCK.FLOAT....",
                                  "X10.AVERAGE_CASP_LUMINESCENCE_UNITS.INT...."),#, 
                                  #"X11.AVERAGE_DAPI_CELL_COUNT.INT....",
                                  #"X12.AVERAGE_DAPI_FIELD_COUNT.FLOAT...."),  
             "replicateCol" = "",
             "plateCol" = "X14.LIBRARY_PLATE.INT....",
             "expCol" = "",
             "ctrlCol" = "",
             "posCtrls" = "",
             "negCtrls" = "",
             "rowCol" = "",
             "colCol" = "",
             "hasCtrls" = "FALSE"
    ))
  }
  else if(datasetName() == "HCC_vorinostat_siRNA"){
    return(list("screenType" = "siRNA",
                "sampleCol" = "X2.ENTREZ_GENE_ID.STRING....gene.target.id", 
                "posColType" = "alpha",
                "posCol" = "X16.WELL_LOCATION.STRING....",
                "accColType" = "Entrez",
                "accCol" = "X2.ENTREZ_GENE_ID.STRING....gene.target.id",
                "measurementCol" = c(#"X5.AVERAGE_CTF_ROBUST_Z_SCORE.FLOAT....", 
                                     #"X6.AVERAGE_CTF_NORMALISED_TO_MOCK.FLOAT....", 
                                     "X7.AVERAGE_CTF_FLUORESCENCE_UNITS.INT....",
                                     #"X8.AVERAGE_CASP_ROBUST_Z_SCORE.FLOAT....",
                                     #"X9.AVERAGE_CASP_NORMALISED_TO_MOCK.FLOAT....",
                                     "X10.AVERAGE_CASP_LUMINESCENCE_UNITS.INT...."),#, 
                                     #"X11.AVERAGE_DAPI_CELL_COUNT.INT....",
                                     #"X12.AVERAGE_DAPI_FIELD_COUNT.FLOAT...."),                                  
                "replicateCol" = "",
                "plateCol" = "X15.LIBRARY_PLATE.INT....",
                "expCol" = "",
                "ctrlCol" = "",
                "posCtrls" = "",
                "negCtrls" = "",
                "rowCol" = "",
                "colCol" = "",
                "hasCtrls" = "FALSE"
    ))
  } 
  else if(datasetName() == "BCSC"){
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
  else if(datasetName() == "KRAS_synleth_compound"){
    return(list("screenType" = "compound",
             "sampleCol" = "Virtual_ID", 
             "posColType" = "alpha",
             "posCol" = "Well",
             "accColType" = "Chembank",
             "accCol" = "ChembankId",
             "measurementCol" = c("RawValueA", "RawValueB", "RawValueC"),
             "replicateCol" = "",
             "plateCol" = "Plate",
             "expCol" = "AssayName",
             "ctrlCol" = "",
             "posCtrls" = "",
             "negCtrls" = "",
             "rowCol" = "",
             "colCol" = "",
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
  else if(datasetName() == "SAVANAH"){
    return(c("screenType" = "",
             "sampleCol" = "Sample", 
             "posColType" = "rowcol",
             "posCol" = "",
             "accColType" = "",
             "accCol" = "Accession",
             "measurementCol" = "PlateReadout",
             "replicateCol" = "Replicate",
             "plateCol" = "Plate",
             "expCol" = "AssayName",
             "ctrlCol" = "",
             "posCtrls" = "",
             "negCtrls" = "",
             "rowCol" = "PlateRow",
             "colCol" = "PlateCol",
             "hasCtrls" = "FALSE"
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
