plot.postprob.hlme <-
function(x,...){
if (!inherits(x, "hlme")) stop("use only with \"hlme\" objects")
if(x$ng==1){
cat("plot.postprob can only be used when ng > 1  \n")
}else{
	if(x$conv==1|x$conv==2) {
		for(i in 1 : x$ng){
			dev.new()
			xlab1 <- paste("class",i)
			title.hist <- paste("distribution of posterior probabilities in class",i)
			hist(x$pprob[,i+2],prob=TRUE,xlab=xlab1,main=title.hist,...)
		}
	}else{
		cat("Output can not be produced since the program stopped abnormally.\n")
	}
}
}



plot.postprob <- function(x,...) UseMethod("plot.postprob")
