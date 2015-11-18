find.mimat <- function(namesToParse){

#we only want to look at mature ids
aliases.mimat <- mirna.aliases[grep("MIMAT", mirna.aliases[,1]),]

#escape * and be sure to catch both mir and miR
namesToParse.pattern <- str_replace(str_replace(as.character(namesToParse), "\\*", "\\\\*"), "mir|miR", "mi(r|R)")

#add semicolon to make sure we don't pick up partial names
namesToParse.pattern <- paste(namesToParse.pattern,";", sep="")

#factors slow things down...
aliases.mimat[,1] <- as.character(aliases.mimat[,1])
aliases.mimat[,2] <- as.character(aliases.mimat[,2])

#search aliases with grep
result <- sapply(namesToParse.pattern, grep, aliases.mimat[,2])

#multiple matches result in NA since we cannot say which one it is
result <- lapply(result, function(x){
  if(length(x) > 1) return(NA)
  else return(aliases.mimat[x[1], 1])
})
return(unlist(as.factor(as.character(result))))
}