observeEvent(input$dataset,{
  
  if(is.null(input$dataset)) return(NULL)
  else if(input$dataset == "none selected") return(NULL)
  else if(input$dataset == "CUSTOM") return(NULL)
  else if(input$dataset == "IMPORTED") return(NULL)
    
  message <- switch(input$dataset,
                    TNFa_Casp4 = "<a href='http://doi.org/10.1128/MCB.06739-11' target='_blank'>Nickles, D., Falschlehner, C., Metzig, M., & Boutros, M. (2012). A Genome-Wide RNA Interference Screen Identifies Caspase 4 as a Factor Required for Tumor Necrosis Factor Alpha Signaling. Molecular and Cellular Biology, 32(17), 3372–3381.</a>",
                    HCC_vorinostat_miRNA = "<a href='http://doi.org/10.1038/sdata.2014.17' target='_blank'>Falkenberg, K. J., Gould, C. M., Johnstone, R. W., & Simpson, K. J. (2014). Genome-wide functional genomic and transcriptomic analyses for genes regulating sensitivity to vorinostat. Scientific Data, 1, 1–13.</a>",
                    HCC_vorinostat_siRNA = "<a href='http://doi.org/10.1038/sdata.2014.17' target='_blank'>Falkenberg, K. J., Gould, C. M., Johnstone, R. W., & Simpson, K. J. (2014). Genome-wide functional genomic and transcriptomic analyses for genes regulating sensitivity to vorinostat. Scientific Data, 1, 1–13.</a>",
                    KRAS_synleth_compound = "<a href='http://chembank.broadinstitute.org/assays/view-project.htm;jsessionid=91940BD33888CDD76E13048EC176BD0D?id=1000525' target='_blank'>ChemBank assay: KRAS Synthetic Lethal Screen in HKE3</a>",
                    "No additional info has been added for this demo dataset")
  showshinyalert(session, "demoDatasetInfo", paste("Reference for the selected demo dataset:<br/><br/> ", message, sep=""), "info")
})

plateMeanInfoFileName <- "help/plateMeanInfo.html"
plateMeanInfoText <- readChar(plateMeanInfoFileName, file.info(plateMeanInfoFileName)$size)

rowColumnInfoFileName <- "help/rowColumnInfo.html"
rowColumnInfoText <- readChar(rowColumnInfoFileName, file.info(rowColumnInfoFileName)$size)

controlPlotInfoFileName <- "help/controlPlotInfo.html"
controlPlotInfoText <- readChar(controlPlotInfoFileName, file.info(controlPlotInfoFileName)$size)

controlPerformancePlotInfoFileName <- "help/controlPerformancePlotInfo.html"
controlPerformancePlotInfoText <- readChar(controlPerformancePlotInfoFileName, file.info(controlPerformancePlotInfoFileName)$size)

zScoreControlPerformancePlotInfoFileName <- "help/zScoreControlPerformancePlotInfo.html"
zScoreControlPerformancePlotInfoText <- readChar(zScoreControlPerformancePlotInfoFileName, file.info(zScoreControlPerformancePlotInfoFileName)$size)

replicateCorrInfoFileName <- "help/replicateCorrInfo.html"
replicateCorrInfoText <- readChar(replicateCorrInfoFileName, file.info(replicateCorrInfoFileName)$size)

normDataTableInfoFileName <- "help/normDataTableInfo.html"
normDataTableInfoText <- readChar(normDataTableInfoFileName, file.info(normDataTableInfoFileName)$size)

wholeScreenInfoFileName <- "help/wholeScreenInfo.html"
wholeScreenInfoText <- readChar(wholeScreenInfoFileName, file.info(wholeScreenInfoFileName)$size)

signalDistributionInfoFileName <- "help/signalDistributionInfo.html"
signalDistributionInfoText <- readChar(signalDistributionInfoFileName, file.info(signalDistributionInfoFileName)$size)

qqPlotInfoFileName <- "help/qqPlotInfo.html"
qqPlotInfoText <- readChar(qqPlotInfoFileName, file.info(qqPlotInfoFileName)$size)

plateViewerInfoFileName <- "help/plateViewerInfo.html"
plateViewerInfoText <- readChar(plateViewerInfoFileName, file.info(plateViewerInfoFileName)$size)

hitsTableInfoFileName <- "help/hitsTableInfo.html"
hitsTableInfoText <- readChar(hitsTableInfoFileName, file.info(hitsTableInfoFileName)$size)

hitsPlotInfoFileName <- "help/hitsPlotInfo.html"
hitsPlotInfoText <- readChar(hitsPlotInfoFileName, file.info(hitsPlotInfoFileName)$size)

heatmapInfoFileName <- "help/heatmapInfo.html"
heatmapInfoText <- readChar(heatmapInfoFileName, file.info(heatmapInfoFileName)$size)

miRNAtargetInfoFileName <- "help/miRNAtargetInfo.html"
miRNAtargetInfoText <- readChar(miRNAtargetInfoFileName, file.info(miRNAtargetInfoFileName)$size)

miRNAhighConfInfoFileName <- "help/miRNAhighConfInfo.html"
miRNAhighConfInfoText <- readChar(miRNAhighConfInfoFileName, file.info(miRNAhighConfInfoFileName)$size)

miRNAfamilyInfoFileName <- "help/miRNAfamilyInfo.html"
miRNAfamilyInfoText <- readChar(miRNAfamilyInfoFileName, file.info(miRNAfamilyInfoFileName)$size)