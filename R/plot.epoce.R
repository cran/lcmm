plot.epoce <- function(x,...){
	if (!inherits(x, "epoce")) stop("use only with \"epoce\" objects")
        dev.new()
	if(x$new.data==FALSE) {
		plot(x$EPOCE[,5]~x$EPOCE[,1],pch=18,lty=1,type="o",ylab="CVPOL",xlab="prediction time",main="Cross-validated prognostic observed log-likelihood",...)
	}else{
		plot(x$EPOCE[,4]~x$EPOCE[,1],pch=18,lty=1,type="o",ylab="MPOL",xlab="prediction time",main="Mean prognostic observed log-likelihood",...)
	}
}

