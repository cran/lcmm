plot.postprob.lcmm <-
function(x,...){
if (!inherits(x, "lcmm")) stop("use only with \"lcmm\" objects")
if(x$conv==1|x$conv==2) {
if(x$ng>1){
for(i in 1 : x$ng){
dev.new()
xlab1 <- paste("class",i)
title.hist <- paste("distribution of posterior probabilities in class",i)
hist(x$pprob[,i+2],prob=TRUE,xlab=xlab1,main=title.hist,...)}}
else  stop("plot only for ng > 1")
}else{
cat("Output can not be produced since the program stopped abnormally.")
stop("Pease check the data & model specification")
}
}



plot.postprob <- function(x,...) UseMethod("plot.postprob")
