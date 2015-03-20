source("source/call_KPM.R")
load("debug/ind.matrix.RData")

ATTACHED_TO_ID <- paste(sample(c(LETTERS[1:6],0:9),32,replace=TRUE),collapse="")

library(foreach)
calls <- foreach(i=1:10) %do%
{
  call.KPM(list(ind.matrix), ATTACHED_TO_ID, Lmin=8, Kmin=1, graphID=2, range=FALSE, with.perturbation=FALSE)  
}
