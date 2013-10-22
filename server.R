require(shiny)
require(plyr)
require(rCharts)
require(gplots)
require(ggplot2)
require(grid)
require(gridExtra)
require(scales)
require(reshape2)
require(RmiR)
require(stringr)
require(VennDiagram)
require(qgraph)

shinyServer(function(input, output) {
  
  source("heatmap.R")
  source("RmiR.R")
  source("go.analysis.R")
  
  #load data
  load.data <- reactive({
    
    if(is.null(input$file)){
      if(input$dataset == "MCF7wt miRNA inhibitor screen")
        data <- read.table("data/MCF7wt.csv", header=T, sep="\t")
      else if(input$dataset == "MCF12Awt miRNA inhibitor screen")
        data <- read.table("data/MCF12Awt.csv", header=T, sep="\t")
    }
    else{
      file <- input$file
      data <- read.table(file$datapath, header=T, sep="\t")
    }
    
    data$Plate <- as.factor(data$Plate)
    data$Replicate <- as.factor(data$Replicate)
    #shorten miRNA label names
    data$Sample <- gsub(",.*", "", gsub("hsa-", "", data$Sample))
    
    return(data)
  })
  
  #filter and summarise
  data <- reactive({

    data <- load.data()

    signalType <- input$normalization
    data$signal <- data[[signalType]]
    input$updateExclusion
    
    data <- isolate({
      if(input$updateExclusion != 0 && nchar(input$exclude) > 0){
        data <- data[-grep(input$exclude, data$miRBase.ID.miRPlus.ID),]
      }
      else data
    })
    deviations <- ddply(data, .(Plate, Well.position), numcolwise(function(x){sd(x)/length(x)}))
    data <- ddply(data, .(Plate, Well.position,miRBase.ID.miRPlus.ID, miRBase.accession, Sample, Control, Comment.in.miRBase, Probe.ID), numcolwise(mean))
    data <- merge(data, deviations[,setdiff(colnames(deviations),c("wellCount", "row", "column", "Probe.ID", "Comment.in.miRBase", "Sample", "Control"))], by=c("Plate", "Well.position"), suffixes=c("", ".sem"))
    return(data)
  })
  
  # find screening hits, e.g. the outliers #
  outliers <- reactive({
    
    input$updateNormalization
    input$updateInclusion
    input$updateExclusion
    
    outl <- isolate({
      outl <- my.outliers(data(), input$method, input$margin, withControls=input$includeControls)
    
    #outl <- ddply(outl, .(Plate, Well.position,miRBase.ID.miRPlus.ID, miRBase.accession, Sample, Control, Comment.in.miRBase, Probe.ID), numcolwise(mean))
    
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
  
  # consensus hit list with all hits
  outliers.all <- reactive({
    hits <- data.frame()
    my.data <- data()
    input$updateNormalization
    
    for(method in input$multiNormalizations)
    {
      my.data$signal <- my.data[[method]]
      outl <- isolate(my.outliers(my.data, input$method, input$margin, withControls=input$includeControls))
      outl$method <- method
      
      outl[outl[[method]] > mean(my.data[[method]], na.rm=T),"category"] <- "promotor"
      outl[outl[[method]] < mean(my.data[[method]], na.rm=T),"category"] <- "suppressor"
      
      hits <- rbind(hits, outl)
    }

    return(ddply(hits, .(Sample, method, category), summarise, count=1))
  })

  # reformat consensus hit list and apply threshold #
  consensusHitList <- reactive({
    consensus <- ddply(outliers.all(), .(Sample, category), summarise, method=paste(method, collapse="/"), count=length(count))
    consensus <- subset(consensus, count >= input$multiThreshold)
    consensus <- merge(consensus, data(), by="Sample")
    consensus <- subset(consensus, select =c(2, 1, seq(3,(ncol(consensus)))))
    consensus$signal <- NULL
    consensus$signal.sem <- NULL
    return(consensus)
  })
  
  ## find mRNA targets ##
  targetsForInteractionGraph <- reactive({
    input$updateTargets
    
    if(input$useConsensus == "hit list") data <- outliers()
    else data <- consensusHitList()
    
    result <- isolate({
     result <- data.frame()
     grouped.targets <- getTargets(data, hits.min=input$at.least.hits, at.least=input$at.least, group.miRNAs=T, group.miRNAs.threshold=input$group.miRNAs.threshold, databases=input$selectedTargetDBs)  
      
     for(i in 1:nrow(grouped.targets)){
        temp.result <- data.frame(
          gene_symbol <- gsub("<|>", "", str_extract(grouped.targets[i, "gene_symbol"], ">.*?<")), 
          miRNA <- str_split(grouped.targets[i, "miRNA_list"], "/")
        )
        
        colnames(temp.result) <- c("gene_symbol", "miRNA")
        result <- rbind(result, temp.result)
      }
      data$miRBase.ID.miRPlus.ID <- sub("mir", "miR", data$miRBase.ID.miRPlus.ID)
      result <- merge(result, unique(data[,c("miRBase.ID.miRPlus.ID", "category")]), by.x="miRNA", by.y="miRBase.ID.miRPlus.ID", all.x=T, all.y=F)
      result[,1] <- gsub("hsa-", "", result[,1])
      print(result)
      return(result)
    })
    
    return(result)
  })
  
  # formatted targets #
  targets <- reactive({
    
    input$updateTargets
    
    if(input$useConsensus == "hit list") data <- outliers()
    else data <- consensusHitList()
    
    result <- isolate({
      result <- getTargets(data, hits.min=input$at.least.hits, at.least=input$at.least, group.miRNAs=input$group.miRNAs, group.miRNAs.threshold=input$group.miRNAs.threshold, databases=input$selectedTargetDBs)
      if(input$excludeDBcol) result <- result[,setdiff(colnames(result), c("db_list"))]
    
      if(nrow(result) == 0) return (data.frame(error="No target has been found. reduce stringency or increase number of hits."))
      if(input$group.miRNAs)
      {
        data$miRBase.ID.miRPlus.ID <- sub("mir", "miR", data$miRBase.ID.miRPlus.ID)
        result$miRNA_list <- as.character(lapply(lapply(str_split(result$miRNA_list, "/"), function(x){
          temp.result <- ""
          for(miR in x){
            if(input$colorizeInTargetList){
              cat.color <- switch(unique(subset(data, miRBase.ID.miRPlus.ID == miR)$category),
                 "included" = "yellow",
                 "suppressor" = "red",
                 "promotor" = "blue") 
              temp.result <- paste(temp.result, "<p style='color:", cat.color, "'>", miR, "</p>", sep="")
            }
            else{
              if(nchar(temp.result) == 0) temp.result <- miR
              else temp.result <- paste(temp.result, "/", miR, sep="")
            } 
          }
          return(temp.result)
        }), paste, collapse=""))
        
        #count promotors (blue) and suppressors (red), build difference
        result$promotor_vs_suppressor <- sapply(result$miRNA_list, function(x){
            blue <- gregexpr("blue", x)[[1]]
            red <- gregexpr("red", x)[[1]]
            if(blue[1] != -1 && red[1] != -1) return(length(blue) - length(red))
            else if(blue[1] == -1) return(-length(red))
            else return(length(blue))
          })
      }
      return(result)
    })                 
    return(result)
  })
  
  ## mirCancerDB ##
  # load database #
  mircancer.database <- reactive({
    read.table("data/miRCancerJune2013.txt", sep="\t", header=T, quote="\"")
  })
  
  # find entries for currently selected hits #
  outliers.mircancer <- reactive({
    outliers <- outliers()
    mirdb <- mircancer.database()
    outliers$miRBase.ID.miRPlus.ID <- sub("hsa-miR", "hsa-mir", outliers$miRBase.ID.miRPlus.ID)
    outliers <- outliers[-grep(".*>NA</a>", outliers$miRBase.accession),]
    merge(outliers[,c("miRBase.ID.miRPlus.ID", "miRBase.accession")], mirdb, by.x="miRBase.ID.miRPlus.ID", by.y="mirId", all.x=T)
  })
  
  ## gene ontology enrichment analysis ##
  goEnrichment <- reactive({
    goEnrichmentAnalysis(targets(), goUpOrDown=input$goUpOrDown, goDomain=input$goDomain, goScoringThreshold=input$goSignThreshold, orderMethod=input$goOrderMethod, topNodes=input$goTopNodes)
  })
  
  ### Plots ###
  
  # plot a miRNA -> mRNA(genes) interaction graph
  output$interactionGraph <- renderPlot({
    
    targets <- targetsForInteractionGraph()
    
    # variable that sets the least amount of miRNA's that it should connect.
    nrConnects <- 2
    
    # Using hashes for the lists, only requires average O(1) time for read/write
    
    #List of genes
    genes <- new.env( hash = TRUE, parent = emptyenv())
    
    #List of interactions including the gene
    geneInteraction <- new.env( hash = TRUE, parent = emptyenv())
    
    #List of category for miRNA
    miRNA_categories <- new.env( hash = TRUE, parent = emptyenv())
    
    print("-----------------")
    
    for(i in 1:nrow(targets)){
      gene_symbol <- as.character(targets[i, "gene_symbol"])
      miRNA <- as.character(targets[i, "miRNA"]) 
      
      # building hash for the number of observed connections, 
      # and the list of genes that it interacted with 
      if(is.null(genes[[gene_symbol]])){
        genes[[gene_symbol]] <- 1
        geneInteraction[[gene_symbol]] <- c(miRNA)
      }else{
        genes[[gene_symbol]] <- genes[[gene_symbol]] + 1     
        geneInteraction[[gene_symbol]] <- append(geneInteraction[[gene_symbol]] , miRNA)
      }
      
      # add the category to list of miRNA categories
      if(is.null(miRNA_categories[[miRNA]])){
        miRNA_categories[[miRNA]] <- targets[i, "category"]
      }
      
      NULL
    }
    
    # The edge list
    mat <- matrix(ncol = 2)
    
    
    # Build the edge list (matrix)
    for(v in ls(genes)){
      # for each gene found, see if the amount of connects is at least the minimum amount
      if(genes[[v]] >= nrConnects){
        for(gI in geneInteraction[[v]]){
          mat <- rbind(mat, c(v, gI))
        }
        
      }
    }
    
    # Must have at least 3 rows (+ 1 for [NA, NA] in matrix), else it fails
    if(nrow(mat) >= 4){ 
      colnames(mat) <- c("from", "to")
      mat = mat[-1,]
      print(mat)
      g <- qgraph(mat,directed=FALSE, layout = "spring")
      plot(g)
    }else{
      plot(1, type="n", axes=F, xlab="", ylab="")
      title(unlist(c("Too few genes with at least ", nrConnects, " connections were found.", "No graph will be drawn.")))
    }
    
  })
  
  # bar chart to visualize how many normalizations report a given hit #
  output$normcomparison <- renderChart({
    data <- outliers.all()
    data <- subset(data, Sample %in% consensusHitList()$Sample)
    p3 <- rPlot(x="Sample", y="count", color='method', data = data, type = 'bar')
    p3$addParams(dom='normcomparison', width="800", height="600")
    return(p3)
  })
  
  # Venn diagram to show overlap of hits in normalization methods #
  output$consensusVennDiagram <- renderPlot({
    data <- outliers.all()
    data <- split(data, data$method)
    if(length(data) <= 5){
      data <- lapply(data, function(x){na.omit(x$Sample)})
      fill.colors <- c("purple", "cornflowerblue", "green", "yellow", "orange", "darkorchid1")
      fill.colors <- fill.colors[1:length(data)]
      p <- venn.diagram(data, filename = NULL, force.unique=T, col="transparent", 
                        fill = fill.colors,
                        alpha=0.5, fontfamily="serif")
      grid.draw(p)
    }
    else textplot("this feature only supports up to 5 normalization methods")
  })
  
  # Heatmap #
  output$heatmapPlot <- renderPlot({
    input$updateNormalization
    isolate(my.heatmap(data(), input$normalization, input$margin, input$method, outliers(), colorA=input$colorA, colorB=input$colorB))
  })
  
  # Interactive plots for screen hits #
  output$scatterPlotHits <- renderChart({
    outl <- outliers()
    p1 <- dPlot(y=input$normalization, x=c("Sample", "Well.position"), z=paste(input$normalization, ".sem", sep=""), data=outl, type="bubble", groups="category", width=1000, height=600)  
    p1$addParams(dom='scatterPlotHits')
    if(input$showSEM) p1$zAxis(type="addMeasureAxis", overrideMax = 2*max(outl[[paste(input$normalization, ".sem", sep="")]]))

    return(p1)
  })
  
  # Plate viewer interactive scatter plot #
  output$scatterPlotI <- renderChart({
     data <- subset(data(), Plate==input$plateSelected)
     p2 <- dPlot(y=input$normalization, x=c("Sample", "Well.position"), z=paste(input$normalization, ".sem", sep=""), data=data, type="bubble", groups="Control", width=1000, height=600)
     p2$xAxis(type="addCategoryAxis", orderRule="wellCount")
     #p2$yAxis(type="addMeasureAxis")
     p2$addParams(dom='scatterPlotI')
     if(input$showSEM) p2$zAxis(type="addMeasureAxis", overrideMax = 2*max(data[[paste(input$normalization, ".sem", sep="")]]))
     return(p2)
  })
  
  # Whole screen scatter plot #
  output$scatterPlot <- renderPlot({
    data <- load.data()
    q <- ggplot(data, aes_string(x="wellCount", y=input$normalization, color="Plate", shape="Replicate")) + geom_point() + geom_line(stat="hline", yintercept="mean", color="black", aes(group=interaction(Plate, Replicate)))
    print(q)
  })
  
  # Plot for control performance (box plots) #
  output$controlPlot <- renderPlot({
    data <- load.data()
    q <- ggplot(subset(data, !is.na(Control)), aes_string(x="Sample", y=input$normalization, color="Plate", shape="Replicate")) + geom_boxplot()
    
    print(q)
  })
  
  # Plot for row and column effects over whole screen #
  output$rowAndColumn <- renderPlot({
    p1 <- qplot(x=row, y=centered, data=load.data(), geom="bar", stat="summary", fun.y="mean", fill=Replicate, position="dodge", main="Row Mean") + scale_y_continuous(labels=percent)
    p2 <- qplot(x=column, y=centered, data=load.data(), geom="bar", stat="summary", fun.y="mean", fill=Replicate, position="dodge", main="Column Mean") + scale_y_continuous(labels=percent)
    print(grid.arrange(p1,p2))
  })
  
  # Plot for comparison of normalization methods in plate viewer #
  output$normalizationComparison <- renderPlot({
    data <- load.data()
    data$signal <- NULL
    p2 <- qplot(x=wellCount, y=value, data=melt(subset(data, Plate==input$plateSelected), id.vars=c("Plate", "Well.position", "wellCount", "Replicate", "Product..", "miRBase.ID.miRPlus.ID", "miRBase.accession", "Sample", "Control", "row", "column", "Probe.ID", "Comment.in.miRBase")), color=Replicate, main="Comparison of different normalization methods") + facet_wrap(~variable, scales="free", ncol=3) + geom_smooth(method="loess")
    print(p2)
  })
  
  #### Generate HTML data tables ###
  
  output$interactionTable <- renderChart2({
    dTable(targetsForInteractionGraph())
  })
  
  # Hit List #
  output$table <- renderChart2({
    data <- outliers()
    data[data$category %in% c("promotor"),"category"] <- "<div style='background:#80B1D3; text-align:center; border-radius: 15px; width:25px; height:25px;'>P</div>"
    data[data$category %in% c("suppressor"),"category"] <- "<div style='background:#FB8072; text-align:center; border-radius: 15px; width:25px; height:25px;'>S</div>"
    data[data$category %in% c("included"),"category"] <- "<div style='background:#FDB462; text-align:center; border-radius: 15px; width:25px; height:25px;'>I</div>"
    dTable(data, sPaginationType='full_numbers')
  })
  
  # Consensus hit list #
  output$consensusHitList <- renderChart2({
    data <- consensusHitList()
    data[data$category %in% c("promotor"),"category"] <- "<div style='background:#80B1D3; text-align:center; border-radius: 15px; width:25px; height:25px;'>P</div>"
    data[data$category %in% c("suppressor"),"category"] <- "<div style='background:#FB8072; text-align:center; border-radius: 15px; width:25px; height:25px;'>S</div>"
    
    dTable(data, sPaginationType='full_numbers')
  })
  
  # mRNA targets #
  output$targets <- renderChart2({
    dTable(targets(), sPaginationType='full_numbers')
  })
  
  # GO enrichment analysis based on mRNA targets #
  output$goEnrichmentTable <- renderChart2({
    if(input$group.miRNAs && input$colorizeInTargetList){
      dTable(goEnrichment(), sPaginationType='full_numbers')
    } else{
      dTable(data.frame(error="This feature only works when miRNAs are grouped and colorized under 'Target Genes'!"))
    }
  })
  
  # mirCancerDB #
  output$mircancerTable <- renderChart2({
    dTable(outliers.mircancer(), sPaginationType='full_numbers')
  })
  
  ### Generate file downloads ###
  
  # Hit list #
  output$downloadHits <- downloadHandler(
    filename = function() { paste('hits', input$normalization, input$margin, input$method, '.csv', sep='_') },
    content = function(file) {
      data <- outliers()
      data$miRBase.url <- gsub("'", "", str_extract(data$miRBase.accession, "'http://.*?'"))
      data$miRBase.accession <- gsub("<|>", "", str_extract(data$miRBase.accession, ">.*?<"))
      data$signal <- NULL
      data$signal.sem <- NULL
      write.table(data, file, row.names=F, sep=",", quote=F)
    }
  )
  
  # Consensus hit list #
  output$downloadConsensusHits <- downloadHandler(
    filename = function() { paste('consensus', 'hits', input$margin, input$method, '.csv', sep='_') },
    content = function(file) {
      data <- consensusHitList()
      data$miRBase.url <- gsub("'", "", str_extract(data$miRBase.accession, "'http://.*?'"))
      data$miRBase.accession <- gsub("<|>", "", str_extract(data$miRBase.accession, ">.*?<"))
      data$signal <- NULL
      data$signal.sem <- NULL
      write.table(data, file, row.names=F, sep=",", quote=F)
    }
  )
  
  # mRNA targets #
  output$downloadTargets <- downloadHandler(
    filename = function() { paste('targets', paste(input$selectedTargetDBs, collapse="_"), input$margin, input$method, '.csv', sep='_') },
    content = function(file) {
      data <- targets()
      if(input$colorizeInTargetList){
        data$categories <- gsub("blue", "P", gsub("red", "S", sapply(str_extract_all(data$miRNA_list, "red|blue"), paste, collapse="/")))
        data$miRNA_list <- gsub("<|>", "", sapply(str_extract_all(data$miRNA_list, ">.*?<"), paste, collapse="/"))
      } 
      data$url <- str_extract(data$gene_symbol, "http://.*?%5D")
      data$gene_symbol <- gsub("<|>", "", str_extract(data$gene_symbol, ">.*?<"))
      
      write.table(data, file, row.names=F, sep=",", quote=F)
    }  
  )
  
  # GO enrichment graph #
  output$dlGnOntGraph <- downloadHandler(
    filename = function() { return("GO.pdf")},
    content = function(file) {
      data <- goEnrichment()
      pdf(file)
      showSigOfNodes(attr(data, "topGO"), score(attr(data, input$goSelectedMethod)), firstSigNodes = input$goSelectedNodes, useInfo=input$goUseInfo)
      dev.off()
    }
  )
  
  # GO enrichment table #
  output$dlGnOntTbl <- downloadHandler(
    filename = function() { return("GO.csv")},
    content = function(file) {
      data <- goEnrichment()
      write.table(data, file, row.names=F, sep=",", quote=F)
    }
  )
  
  output$datasetName <- renderPrint({
    input$updateNormalization
    isolate(input$dataset)
  })
}) 
