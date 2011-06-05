

plot.baselinerisk.Jointlcmm <- function(x,legend.loc="topleft",ylim=NULL,...){

if(missing(x)) stop("The argument x should be specified")
if (!inherits(x, "Jointlcmm")) stop("use only with \"Jointlcmm\" objects")

if((x$conv==1|x$conv==2)& (sum(is.na(x$predSurv)==0))){
title1 <- "Class-specific baseline risk functions"

ylim1 <- ylim
if (is.null(ylim)){
ylim1 <- c(min(x$predSurv[,2:x$ng]),max(x$predSurv[,2:x$ng]))
}

for (g in 1:x$ng){
Y <- x$predSurv[,(1+g)]
plot(Y~x$predSurv[,1],col=g,type="l",ylim=ylim1,xlim=c(min(x$predSurv[,1]),max(x$predSurv[,1])),xlab="Time",pch=2,ylab="Baseline risk function",main=title1,lty=g)
par(new=TRUE)
}

temp <- paste("class",1:x$ng,sep="") 
legend(x=legend.loc,legend=temp,col=1:x$ng,lty=1:x$ng,inset=.02)
par(new=FALSE)
}else{
cat("Output can not be produced. The program stopped abnormally or there was an error in the computation of the estimated baseline risk functions and survival functions.\n")
}
}


plot.baselinerisk <- function(x,legend.loc="topleft",ylim=NULL,...) UseMethod("plot.baselinerisk")

