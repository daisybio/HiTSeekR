source("source/call_KPM.R")
load("debug/ind.matrix.RData")

ATTACHED_TO_ID <- paste(sample(c(LETTERS[1:6],0:9),32,replace=TRUE),collapse="")
url <- "http://tomcat.compbio.sdu.dk/keypathwayminer/"

library(foreach)
calls <- foreach(L=1:1) %do%
{
  call.KPM(list(ind.matrix), ATTACHED_TO_ID, Lmin=10, Kmin=1, graphID=5, range=FALSE, with.perturbation=FALSE, url=url)  
}
