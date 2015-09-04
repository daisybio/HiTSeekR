#download diana microT CDS dataset after registering at the diana lab website
#for every miRNA-mRNA interaction the dataset contains extra lines indicating where in the sequence 
#the interaction occurs. these lines start with either CDS or UTR3 and can be filtered efficiently with 

#awk -F, '$1 !~ /(^CDS|^UTR3)/ && $3 ~ /^hsa/' microT_CDS_data.csv > microT_CDS_data.filtered.hsa.csv

#in this process, I also filtered for human miRNAs. Next, we split some columns

#read file efficiently
library(sqldf)
microT.hsa <- read.csv.sql(file="data//microT_CDS_data.filtered.hsa.csv",header=T, dbname="microT.hsa.sqlite3", sql="select * from file", row.names=F,field.types=list(TranscriptId="text", "GeneId(name)"="text", "Mirna-Name(miRBase-version)"="text", "miTG-score"="real"), colClasses=c("character", "character", "character", "numeric"))
gene_split <- colsplit(microT.hsa[,2], "\\(", c("ensembl_id", "symbol"))

mirna_split <- colsplit(microT.hsa[,3], "\\(", c("mirna_id", "miRBase_version"))
score <- microT.hsa[,"miTG_score"]

microT.hsa.processed <- data.frame(gene_split, mirna_split, score)

#remove extra ) brackets
microT.hsa.processed <- within(microT.hsa.processed, symbol <- substr(symbol, 1, nchar(symbol)-1))
microT.hsa.processed <- within(microT.hsa.processed, miRBase_version <- as.integer(substr(miRBase_version, 1, nchar(miRBase_version)-1)))

#add some more info from mirbase
library(mirbase.db)
miRNA_acc_table <- as.data.frame(mirbaseID2ACC)
microT.hsa.processed <- left_join(microT.hsa.processed, miRNA_acc_table, by=c("mirna_id"))

miRNA_mat_table <- as.data.frame(mirbaseMATURE)[,-7] #omit experiment type
library(dplyr)
microT.hsa.processed <- left_join(microT.hsa.processed, miRNA_mat_table, by=c("mirna_id"))

#add entrez gene ids
library(org.Hs.eg.db)
gene_acc_table <- as.data.frame(org.Hs.egENSEMBL2EG)
microT.hsa.processed <- left_join(microT.hsa.processed, gene_acc_table, by=c("ensembl_id"))
