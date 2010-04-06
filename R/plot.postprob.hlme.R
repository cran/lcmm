plot.postprob.hlme <-
function(x,...){
if (!inherits(x, "hlme")) stop("use only with \"hlme\" objects")
if(x$ng>1){
for(i in 1 : x$ng){
dev.new()
xlab1 <- paste("class",i)
title.hist <- paste("distribution of posterior probabilities in class",i)
hist(x$pprob[,i+1],prob=TRUE,xlab=xlab1,main=title.hist,...)}}
else  stop("plot only for ng > 1")
}



plot.postprob <- function(x,...) UseMethod("plot.postprob")
