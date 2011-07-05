plot.postprob.Jointlcmm <-
function(x,...){
if (!inherits(x, "Jointlcmm")) stop("use only with \"Jointlcmm\" objects")
if(x$ng==1){
cat("plot.postprob can only be used when ng > 1  \n")
}else{
	if(x$conv==1|x$conv==2) {
		if(sum(is.na(x$pprob[,3:(x$ng+2)]))==0){
			for(i in 1 : x$ng){
				dev.new()
				xlab1 <- paste("class",i)
				title.hist <- paste("Posterior probabilities given longitudinal and time-to-event data in class",i)
				hist(x$pprob[,i+2],prob=TRUE,xlab=xlab1,main=title.hist,...)
			}
		}else{
		cat("Error in the computation of posterior class-membership probabilities given all the information \n")
		}	       
   }else{
      cat("Output can not be produced since the program stopped abnormally.\n")
   }
}
}


plot.postprob <- function(x,...) UseMethod("plot.postprob")
