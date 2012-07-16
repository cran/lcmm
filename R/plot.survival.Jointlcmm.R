

plot.survival.Jointlcmm <- function(x,legend.loc="topright",ylim=c(0,1),...){

if(missing(x)) stop("The argument x should be specified")
if (!inherits(x, "Jointlcmm")) stop("use only with \"Jointlcmm\" objects")
if((x$conv %in% c(1,2)) & (sum(is.na(x$predSurv)==0))){

dev.new()
title1 <- "Class-specific event-free probability"
for (g in 1:x$specif[[3]]){
Y <- exp(-x$predSurv[,(1+x$specif[[3]]+g)])
plot(Y~x$predSurv[,1],col=g,type="l",ylim=ylim,xlim=c(min(x$predSurv[,1]),max(x$predSurv[,1])),xlab="Time",pch=2,ylab="Event-free probability",main=title1,lty=g,...)
par(new=TRUE)
}

temp <- paste("class",1:x$specif[[3]],sep="") 
legend(x=legend.loc,legend=temp,col=1:x$specif[[3]],lty=1:x$specif[[3]],inset=.02,bty="n")
par(new=FALSE)
}else{
cat("Output can not be produced. The program stopped abnormally or there was an error in the computation of the estimated baseline risk functions and survival functions.\n")
}
}


plot.survival <- function(x,legend.loc="topright",ylim=c(0,1),...) UseMethod("plot.survival")
