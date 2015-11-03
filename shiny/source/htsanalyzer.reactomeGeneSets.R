library(reactome.db)
library(dplyr)

getReactomeGeneSets <- function(species="Homo sapiens"){
gene.sets <- as.data.frame(reactomeEXTID2PATHID)
pathway.names <- as.data.frame(reactomePATHID2NAME)

#group by pathway ID and collect sets as list
gene.sets.grouped <- gene.sets %>% group_by(DB_ID) %>% dplyr::summarize(gene.sets = list(gene_id))

#add pathway names
named.gene.sets <- left_join(gene.sets.grouped, pathway.names, by="DB_ID")

#filter for specific species via pathway name
if(!is.null(species))
named.gene.sets <- named.gene.sets %>% filter(grepl(species, path_name))

#remove duplicated pathway names for the same ID (should only appear with cross species)
named.gene.sets <- named.gene.sets[!duplicated(named.gene.sets$DB_ID),]
list.of.gene.sets <- named.gene.sets$gene.sets
names(list.of.gene.sets) <- named.gene.sets$DB_ID
attr(list.of.gene.sets, "pathway_names") <- named.gene.sets[,c("DB_ID", "path_name")]
return(list.of.gene.sets)
}