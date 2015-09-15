# #faster standard deviation calculation
# Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
# cppFunction(' double sdC(NumericVector v){
# double sum = std::accumulate(std::begin(v), std::end(v), 0.0);
# double m =  sum / v.size();
# 
# double accum = 0.0;
# std::for_each (std::begin(v), std::end(v), [&](const double d) {
#     accum += (d - m) * (d - m);
# });
# 
# double stdev = sqrt(accum / (v.size()-1));
# return(stdev);              
#               }')  

#filter and summarise
data <- reactive({
  #load data in "normalized" form with known column names
  data <- processedData()
  
  if(is.null(data)) stop("Please process the input data first")
  #fire up a progress bar
  progress <- shiny::Progress$new()
  progress$set(message = "Merge & Filter...", value = 0)
  on.exit(progress$close())
  
  #update progress bar function
  updateProgress <- function(value = NULL, detail = NULL) {
    if (is.null(value)) {
      value <- progress$getValue()
      value <- value + (progress$getMax() - value) / 5
    }
    progress$set(value = value, detail = detail)
  }
  
  #update progress bar
  updateProgress(detail = "Merging replicates", value=0.2)

  merge.funs <- funs(mean(., na.rm=T), sd(., na.rm=T))
  
  #data <- data %>% filter(!is.na(Raw))
  data <- data %>% group_by(Experiment, Readout, Plate, Row, Column, Sample, Accession, Well.position, Control) %>% summarise_each(merge.funs, c(-Replicate, -wellCount))  
  
  #fix column names and get rid of NaNs
  colnames(data) <- str_replace(colnames(data), "_mean", "")
  
  #if we are dealing with miRNAs add miRNA family name and ID
  if(input$screenType == "miRNA")
  {
    #update progress bar
    progress$set(message = "Querying mirbase", value=0.6)
    
    fam <- as.data.frame(mirbaseFAMILY)
    
    mergeRows <- function(y){      
      if(length(unique(y)) > 1) return(paste(y, collapse="/"))
      else return(y[1])
    }     
    data$Accession <- as.character(data$Accession)
    
    if(input$accessionColType=="mature_name")
    {
      updateProgress(detail = "Replace alias with MIMAT", value=0.7)
      data$Accession <- find.mimat(data$Accession)
    }
    if(input$accessionColType %in% c("MIMAT", "mature_name"))
    {
      updateProgress(detail = "Get miRNA family id", value=0.8)
      mimat <- as.data.frame(mirbaseMATURE)
      result <- left_join(mimat, fam, by="mirna_id")
      data <- left_join(data, result, by=c("Accession" = "mature_acc"))
      
      data <- data %>% group_by(Experiment, Readout, Plate, Well.position) %>% summarise_each(funs(mergeRows))      
    }    
    else if(input$accessionColType=="MI"){
      updateProgress(detail = "Get miRNA id", value=0.8)
      mi <- as.data.frame(mirbaseACC2ID)
      data <- left_join(data, mi, by=c("Accession" = "mirna_acc"))
      data <- left_join(data, fam, by=c("mirna_id"))
    }    
    
    else{
      stop("unknown miRNA identifier selected")
    }
  }
  else if(input$screenType == "siRNA")
  {
    if(input$accessionColType=="FlybaseCG")
    {
      library(org.Dm.eg.db)
      flybaseCG <- as.data.frame(org.Dm.egFLYBASECG)
      data <- left_join(data, flybaseCG, by=c("Accession" = "flybase_cg_id"))
    }
    else if(input$accessionColType == "Entrez"){
      data$gene_id <- as.character(data$Accession)
    }
    else if(input$accessionColType == "Ensembl"){
      ensembl <- as.data.frame(org.Hs.egENSEMBL2EG)
      data <- dplyr::left_join(data, ensembl, by=c("Accession" = "ensembl_id")) %>% dplyr::rename(gene_symbol=alias_symbol)      
    }
    else if(input$accessionColType == "GeneSymbol"){
      aliassymbol <- as.data.frame(org.Hs.egSYMBOL)
      data <- dplyr::left_join(data, aliassymbol, by=c("Accession" = "symbol"))       
    }
    
    else{
      stop("accession type not supported")
    }
    if(input$accessionColType != "FlybaseCG"){      
      genesymbol <- as.data.frame(org.Hs.egSYMBOL)
      symbols <- dplyr::ungroup(data) %>% dplyr::select(gene_id)
      symbols <- dplyr::left_join(symbols, genesymbol, by=c("gene_id")) %>% dplyr::rename(gene_symbol=symbol)           
      symbols <- symbols %>% dplyr::distinct() %>% dplyr::group_by(gene_id) %>% dplyr::summarise(gene_symbol = paste(gene_symbol, collapse="/"))
      data <- dplyr::left_join(data, symbols, by="gene_id")
    }
  }
  
  return(as.data.frame(data))
})
