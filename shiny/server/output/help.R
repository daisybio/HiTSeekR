plateMeanInfoFileName <- "help/plateMeanInfo.html"
plateMeanInfoText <- readChar(plateMeanInfoFileName, file.info(plateMeanInfoFileName)$size)

rowColumnInfoFileName <- "help/rowColumnInfo.html"
rowColumnInfoText <- readChar(rowColumnInfoFileName, file.info(rowColumnInfoFileName)$size)

controlPlotInfoFileName <- "help/controlPlotInfo.html"
controlPlotInfoText <- readChar(controlPlotInfoFileName, file.info(controlPlotInfoFileName)$size)

controlPerformancePlotInfoFileName <- "help/controlPerformancePlotInfo.html"
controlPerformancePlotInfoText <- readChar(controlPerformancePlotInfoFileName, file.info(controlPerformancePlotInfoFileName)$size)

replicateCorrInfoFileName <- "help/replicateCorrInfo.html"
replicateCorrInfoText <- readChar(replicateCorrInfoFileName, file.info(replicateCorrInfoFileName)$size)