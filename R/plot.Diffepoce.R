plot.Diffepoce <- function(x,ylim=NULL,xlim=NULL,...){
	if (!inherits(x, "Diffepoce")) stop("use only with \"Diffepoce\" objects")

	if (is.null(ylim)&all(is.na(x$DiffEPOCE[,3]))|all(is.na(x$DiffEPOCE[,4]))) stop("can't produce the plot with missing differences in EPOCE")


      dev.new()

	if(x$new.data==FALSE) {
		title <- expression(Delta(CVPOL))
	}else{
		title <- expression(Delta(MPOL))	
	}

ylim1 <- ylim
if (is.null(ylim)){
ylim1 <- c(min(x$DiffEPOCE[!(is.na(x$DiffEPOCE[,2])),2],x$DiffEPOCE[!(is.na(x$DiffEPOCE[,3])),3]),max(x$DiffEPOCE[!(is.na(x$DiffEPOCE[,2])),2],x$DiffEPOCE[!(is.na(x$DiffEPOCE[,4])),4]))	
}
xlim1 <- xlim
if (is.null(xlim)){
xlim1 <- c(min(x$DiffEPOCE[,1]),max(x$DiffEPOCE[,1]))
}

plot(x$DiffEPOCE[,1],x$DiffEPOCE[,2],pch=18,lty=1,type="o",xlab="prediction time",ylab=title,main="Difference in EPOCE estimates",bty="l",ylim=ylim1,xlim=xlim1,...)
par(new=T)
plot(x$DiffEPOCE[,1],x$DiffEPOCE[,3],pch=18,type="o",lty=3,xlab="prediction time",ylab=title,main="Difference in EPOCE estimates",bty="l",ylim=ylim1,xlim=xlim1,...)
par(new=T)
plot(x$DiffEPOCE[,1],x$DiffEPOCE[,4],pch=18,type="o",lty=3,xlab="prediction time",ylab=title,main="Difference in EPOCE estimates",bty="l",ylim=ylim1,xlim=xlim1,...)
par(new=T)
plot(x$DiffEPOCE[,1],rep(0,length(x$DiffEPOCE[,1])),lty=1,type="l",col="lightgrey",xlab="prediction time",ylab=title,bty="l",ylim=ylim1,xlim=xlim1,...)

}



