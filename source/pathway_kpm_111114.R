# blame for bugs: Richa
# last modified: 09-05-14
library(gtools)
require (ggplot2) # for line plots
library (grid) # for plotting in grids
library (useful) # vplayout

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(plotlist, cols, layout) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- plotlist
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Datasets to be used with KPM
# read the dataset input provided in input parameters
# convert it into command line argument that be used by stand alone KPM along with other parameters
# also writes a dataset file
# percentage of patient cases are computed if specified else it puts l as it is
GetDatasetFile <- function (datasets, l, data.file, perc) {
 datasets <- unlist(strsplit(as.character(datasets), " "))
 num.cases <- unlist(lapply(datasets, FUN = function (x) ncol(read.delim(x)) - 1))
 data.set <- matrix(nrow = length(datasets), ncol = 3)
 data.set <- data.frame(data.set)
 names(data.set) <- c("ID", "L", "fileloc")
 data.set$fileloc <- datasets
 if (perc == 1) {
  data.set$L <- floor(num.cases*(l*0.01))
  small.l <- which(num.cases < 10) # can be an user input
  if(length(small.l) > 0){
   data.set$L[small.l] <- 0
  }
 } else {
    data.set$L <- l
 }
 data.set$ID <- 1:length(datasets)
 mats <- (unlist(sapply(1:nrow(data.set), FUN = function (x) paste(" -matrix", x, "=", data.set$fileloc[x], sep = ""))))
 ls <- (unlist(sapply(1:nrow(data.set), FUN = function (x) paste(" -L", x, "=", data.set$L[x], sep=""))))
 data.arg <- paste(mats, ls, collapse = " ")
 write.table(data.set, file = data.file, quote = F, sep = "\t", col.names = T, row.names = F)
 return(data.arg)
}

# Run KPM over a set of datasets, ls, ks
RunKpm <- function (algorithm, strategy, name.run, datasets, lowL, upL, stepL, kpm.loc, lowK, upK, stepK, net, combine.operation, combine.formula, 
                    results.folder, percentage) {
 name.run <- paste(algorithm, strategy, name.run, sep = "_")
 ## loop over (1) datasets, (2) k and (3) l values
 for (l in seq(from = lowL, to = upL, by = stepL)){
  data.arg <- GetDatasetFile (datasets = datasets, l, data.file = paste(kpm.loc, "datasets_file.txt", sep = ""), percentage)
  print(data.arg)
  for (k in seq(from = lowK, to = upK, by = stepK)){
   kpm.ex <- paste("java -jar ", kpm.loc, "/KPM-4.0.jar", " -program=KPM", " -algo=", algorithm, 
                   " -strategy=", strategy, " -K=", k, " -suffix=", paste(name.run, "_K", k, "_L", l, sep = ""), 
                   " -graphFile=", net, " -combineOp=", combine.operation, " -combineFormula=", "'", combine.formula, "'", 
                   " -resultsDir=", results.folder, data.arg, sep = "")
  system (kpm.ex)
  print (kpm.ex)
  }
 }
}

# given the pathways computed by KPM and the gold standard file 
# it computed a matrix with k, l, ppr, size of the subnetwork
GetKLMat <- function (files, ref.file) {
 # get unique ks and ls from the list of files
 ls <- unique(unlist(lapply(files, FUN = function (x) unlist(strsplit(unlist(strsplit(as.character(x), "_L"))[2], ".txt"))[1])))
 ks <- unique(unlist(lapply(files, FUN = function (x) unlist(strsplit(unlist(strsplit(as.character(x), "_K"))[2], "_"))[1])))
 # create a matrix that will store the k, l and the ratio for all the pathways
 mat <- matrix(nrow = (length(ks) * length(ls)), ncol = 5)
 mat <- data.frame(mat)
 names(mat) <- c("k","l","ratio", "ratio2", "size")
 # for each K and L read the pathways from the list of files
 row.num <- 0
 for (i in c(1:length(ks))){
  for (j in c(1:length(ls))){
   files.ij <- files[grep(paste("K", ks[i], "_", "L", ls[j], sep = ""), files)]
   path.nodes <- unique(unlist(lapply(files.ij, FUN = function (x) {pn <- read.delim(x); unique(c(as.matrix(pn[, 1]), as.matrix(pn[, 2])))})))
   # get k and L value that we added as a sufix during the run
   row.num <- row.num + 1
   mat$k[row.num] <- ks[i]
   mat$l[row.num] <- ls[j]
   if (file.exists(ref.file)) {
    # computing the ratio ref_genes_in_pathway/total_genes_in_pathway
    ref.genes <- read.delim(ref.file)$genes
    mat$ratio[row.num] <- length(which(path.nodes %in% ref.genes)) / length(path.nodes)
    mat$ratio2[row.num] <- length(which(ref.genes %in% path.nodes)) / length(ref.nodes)
   }
   mat$size[row.num] <- length(path.nodes)   
  }
 }
 # order the matrix as per K values
 mat$k <- as.numeric(mat$k)
 mat <- mat[order(mat$k), ]
 # factorise the ls
 mat$l <- factor(mat$l, levels = c(unique(mat$l)), labels = c(sapply(unique(mat$l), FUN = function (x) x)))
 return(mat)
}

# given the output generated by GetKLMat it plots the ppr plots for each l
GetEvaluationRocPerL <- function (mat) {
 ls <- sort(as.numeric(as.matrix(unique(mat$l))))
 ppr.graphs.list <- list()
 size.graphs.list <- list()
 ppr2.graphs.list <- list()
 num.plot = 1
 for (j in c(1:length(ls))){
  print(j)
  mati <- mat[which(mat[, 2]%in%ls[j]), ]
  mati$k <- as.numeric(mati$k)
  mati$size <- as.numeric(mati$size)
  mati <- mati[order(mati$k), ]
  ppr.graphs.list[[num.plot]] <- ggplot(data = mati[, c(1, 3)], aes(x = k, y = ratio)) + geom_point() + labs (title = paste("positive prediction rate at : L-", ls[j], sep = "" ), x = ("Gene exceptions (K)"), 
                                                                                                          y = ("percentage reference nodes in a pathway"))
  ppr2.graphs.list[[num.plot]] <- ggplot(data = mati[, c(1, 4)], aes(x = k, y = ratio2)) + geom_point() + labs (title = paste("positive prediction rate at : L-", ls[j], sep = "" ), x = ("Gene exceptions (K)"), 
                                                                                                              y = ("fraction of reference nodes"))
  size.graphs.list[[num.plot]] <- ggplot(data = mati[, c(1, 5)], aes(x = k, y = size)) + geom_point() + labs (title = paste("network size acceraltion at : L-", ls[j], sep = "" ),
                                                                                                         x = ("Gene exceptions (K)"), y = ("number of nodes"))
  num.plot <- num.plot + 1
 }
 multiplot(ppr.graphs.list, 2 , layout = NULL)
 multiplot(size.graphs.list, 2 , layout = NULL)
 multiplot(ppr2.graphs.list, 2 , layout = NULL)
}

# given the directory where KPM pathways are kept, the gold standard file, top X number of pathways to be combined, a output file name
# plots the ppr rate for each l over k as x axis using GetEvaluationRocPerL
# also plots a single plot with all ls in different color for a global prespective
GetEvaluationRocTopX <- function (filename, ref.file, results.dir, topX) {
  pdf(file = filename, height = 10, width = 10)
  # controls the font size of the plot
  opar <- par(ps = 12)
  # read the ref gene file whose column name with ref genes should be "genes"
  ref.genes <- read.delim(ref.file, row.names=NULL)$genes
  # get top X pathways files
  files <- unlist(lapply(1:topX, FUN = function (x) unlist(list.files(path = results.dir, pattern = paste("Pathway-0", x, "-INTERACTIONS", sep=""), 
                                                                      full.names=T))))
  # remove the files that are empty
  files <- unlist(lapply(files, FUN = function (x) {if (file.info(x)$size > 0){x}}))
  mat <- GetKLMat(files, ref.file)
  GetEvaluationRocPerL (mat)
  grid.newpage ()
  pushViewport (viewport (layout = grid.layout(2, 2)))
  p <- qplot(k, ratio, data = mat, geom = c("point", "line"), formula = y~x, color = l, 
             main = "positive prediction rate", xlab="node exceptions", 
             ylab = "% reference nodes in the top pathway")
  p2 <- qplot(k, ratio2, data = mat, geom = c("point", "line"), formula = y~x, color = l, 
             main = "positive prediction rate2", xlab="node exceptions", 
             ylab = "fraction of reference nodes in the top pathway")
  print (p, vp = vplayout(1, 1))
  print (p2, vp = vplayout(1, 2))
  opar
  dev.off ()
  save (mat, file = paste(filename, "_.RData", sep = ""))
}

# Finds out optimal K for each of the l used
# plots regression for each l
# plots optimal k for each l versus l
KOptPlot <- function (filename, results.dir) {
 files <- unlist(list.files(path = results.dir, pattern = "Pathway-01-INTERACTIONS", full.names=T))
 files <- unlist(lapply(files, FUN = function (x) {if (file.info(x)$size > 0){x}}))
 mat <- GetKLMat(files, ref.file = "")
 ls <- mixedsort((unique(mat$l)))
 graphs.list <- as.list(1:length(ls))
 k.opt <- matrix(nrow = length(ls), ncol = 2)
 k.opt <- data.frame(k.opt)
 names(k.opt) <- c("k", "l")
 k.opt[,2] <- as.matrix(ls )
 
 for (j in c(1:length(ls))){
  print(j)
  mati <- mat[which(mat[, 2]%in%ls[j]), ]
  mati$k <- as.numeric(mati$k)
  mati$size <- as.numeric(mati$size)
  mati <- mati[order(mati$k), ]
  graphs.list[[j]] <- qplot(k, size, data = mati, geom = c("point"), formula = y~x,
              main = paste("subnetwork size distribution at ", ls[j], " percent", sep = ""), xlab = "node exceptions(K)",
               ylab = "number of nodes in biggest pathway") + stat_smooth(method = "lm", se = FALSE)
  fit <- lm(size~k, data = mati)
  k.opt[j,1] <- mati$k[which(fit$residuals == max(fit$residuals))]
 }
 
 pdf(file = filename, height = 10, width = 10)
 opar <- par(ps = 12)
 multiplot(graphs.list, 2 , layout = NULL)
 grid.newpage ()
 pushViewport (viewport (layout = grid.layout(2, 2)))
 p <- qplot(l, k, data = k.opt, geom = c("point", "line"), formula = y~x, 
            main = "", xlab="percentage case exceptions", 
            ylab = "optimal k")
 print (p, vp = vplayout(1, 1))
 opar
 dev.off()
}

# computes composite subnetwork by combining the specified number of top x pathways
GetCompositeKeyPathway <- function (results.dir, topX) {
 files <- unlist(lapply(1:topX, FUN = function (x) unlist(list.files(path = results.dir, pattern = paste("Pathway-0", x, "-INTERACTIONS", sep=""), 
                                                                      full.names=T))))
 # remove the files that are empty
 files <- unlist(lapply(files, FUN = function (x) {if (file.info(x)$size > 0){x}}))
 # get unique ks and ls from the list of files
 ls <- sort(as.numeric(unique(unlist(lapply(files, FUN = function (x) unlist(strsplit(unlist(strsplit(as.character(x), "_L"))[2], "_"))[1])))))
 ks <- sort(as.numeric(unique(unlist(lapply(files, FUN = function (x) unlist(strsplit(unlist(strsplit(as.character(x), "_K"))[2], "_"))[1])))))
 path.edges <- as.list(length(ks)*length(ls))
 k <- 1
 path.names <- NULL
 for (i in c(1:length(ks))) {
  for (j in c(1:length(ls))) {
   files.ij <- files[grep(paste("K", ks[i], "_", "L", ls[j], sep = ""), files)]
   path.edges[[k]] <- unique(do.call(rbind,lapply(files.ij, FUN = function (x) {pn <- read.delim(x); as.matrix(pn) })))
   k <- k+1
   path.names <- c(path.names, paste("K", ks[i], "_", "L", ls[j], sep = ""))
  }
 }  
 names(path.edges) <- path.names
 return(path.edges)
}
GetEvaluationRocTopXReps <- function (filename, ref.file, results.dir, topX) {
  pdf(file = filename, height = 10, width = 10)
  # controls the font size of the plot
  opar <- par(ps = 12)
  # read the ref gene file whose column name with ref genes should be "genes"
  ref.genes <- read.delim(ref.file, row.names=NULL)$genes
  # get top X pathways files
  files <- unlist(lapply(1:topX, FUN = function (x) unlist(list.files(path = results.dir, pattern = paste("Pathway-0", x, "-INTERACTIONS", sep=""), 
                                                                      full.names=T))))
  # remove the files that are empty
  files <- unlist(lapply(files, FUN = function (x) {if (file.info(x)$size > 0){x}}))
  mat <- GetKLMat(files, ref.file)
  grid.newpage ()
  pushViewport (viewport (layout = grid.layout(2, 2)))
  p <- qplot(k, ratio, data = mat, geom = c("point", "line"), formula = y~x, color = l, 
             main = "positive prediction rate", xlab="node exceptions", 
             ylab = "% reference nodes in the top pathway")
  print (p, vp = vplayout(1, 1))
  opar
  dev.off ()
  save (mat, file = paste(filename, "_.RData", sep = ""))
}


CallKPM <- function(kpm.wrapper, kpm.loc, name.run, ref.file, results.folder, algorithms=c("GREEDY"), strategies=c("INES"),
                    datasets,combine.formula = "(L1 || L2)", combine.operation = "OR", percentage=0, lowL=0, upL=0, stepL=0,
                    lowK=0, upK=50, stepK=10, net) {
 source (kpm.wrapper)
 setwd(kpm.loc)
 if(!file.exists(results.folder)) {
  dir.create(path = results.folder, recursive = T)
 } 
 
 # call kpm
 lapply(algorithms, FUN = function (a) lapply(strategies,FUN=function(s)RunKpm(algorithm = a, strategy = s, name.run, datasets, lowL,
      upL, stepL, kpm.loc, lowK, upK, stepK, net, combine.operation, combine.formula, results.folder, percentage)))
 
 # separate files for each run and then run plots for each of these
 all_files <- list.files(path = results.folder, pattern=".txt")
 folders <- unlist(lapply(1:length(algorithms), FUN = function (x) 
   unlist(lapply(1:length(strategies), FUN = function(y)paste(algorithms[x], strategies[y], sep = "_")))))
                                                       
 lapply(folders, FUN = function (x) {
  dir.create(path = paste(results.folder, "/", x, sep = ""));
  filesx = list.files(path = results.folder, pattern = x);
  lapply(filesx, x, FUN = function(y, x) 
             file.rename(from = paste(results.folder, "/", y, sep = ""), to = paste(results.folder, "/", x, "/", y, sep = "")))})
 
 # plotting the evaluation curve
 topXs <- c(1, 3, 5)
 sapply(topXs, function (topX)
  lapply(folders, FUN = function(x)
                          GetEvaluationRocTopX(
                           filename = paste(results.folder, "/TopX_Evaluation_", x, "_", topX, ".pdf", sep = ""),
                           ref.file, results.dir = paste(results.folder, "/", x, sep = ""), topX)
  ))
}
CallPlots <- function(plot.wrapper, ref.file, results.folder, algorithms=c("GREEDY"), strategies=c("INES")) {
  source(plot.wrapper)
  # separate files for each run and then run plots for each of these
  folders <- unlist(lapply(1:length(algorithms), FUN = function (x) 
    unlist(lapply(1:length(strategies), FUN = function(y)paste(algorithms[x], strategies[y], sep = "_")))))
  
  # plotting the evaluation curve
  topXs <- c(1, 3, 5)
  sapply(topXs, function (topX)
    lapply(folders, FUN = function(x)
      GetEvaluationRocTopX(
        filename = paste(results.folder, "/ppr2_TopX_Evaluation_", x, "_", topX, ".pdf", sep = ""),
        ref.file, results.dir = paste(results.folder, "/", x, sep = ""), topX)
    ))
}
