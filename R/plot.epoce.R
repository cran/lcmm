plot.epoce <- function(x,ylim=NULL,xlim=NULL,...){
	if (!inherits(x, "epoce")) stop("use only with \"epoce\" objects")

	if (is.null(ylim)&all(is.na(x$EPOCE[,4]))|all(is.na(x$EPOCE[,5]))) stop("can't produce the plot with missing EPOCE")

        dev.new()
	if(x$new.data==FALSE) {

ylim1 <- ylim
if (is.null(ylim)){
ylim1 <- c(min(x$EPOCE[!(is.na(x$DiffEPOCE[,5])),5],x$EPOCE[!(is.na(x$DiffEPOCE[,4])),4]),max(x$EPOCE[!(is.na(x$EPOCE[,4])),4],x$EPOCE[!(is.na(x$DiffEPOCE[,5])),5]))	
}
xlim1 <- xlim
if (is.null(xlim)){
xlim1 <- c(min(x$EPOCE[,1]),max(x$EPOCE[,1]))
}


		plot(x$EPOCE[,5]~x$EPOCE[,1],pch=18,lty=1,type="o",ylab="CVPOL",xlab="prediction time",main="Cross-validated prognostic observed log-likelihood",...)
	}else{
		plot(x$EPOCE[,4]~x$EPOCE[,1],pch=18,lty=1,type="o",ylab="MPOL",xlab="prediction time",main="Mean prognostic observed log-likelihood",...)
	}
}

