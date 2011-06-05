plot.linkfunction.lcmm <- function(x,legend.loc="topright",...){

if(missing(x)) stop("The argument x should be specified")
if (!inherits(x, "lcmm")) stop("use only with \"lcmm\" objects")
if(x$conv %in% c(1,2)){
title1 <- "Estimated link function"
plot(x$estimlink[,1]~x$estimlink[,2],type="l",ylim=c(min(x$estimlink[,1]),max(x$estimlink[,1])),xlab="Latent process",pch=2,ylab="Longitudinal outcome",main=title1)
par(new=FALSE)
}else{
cat("Output can not be produced since the program stopped abnormally.")
stop("Pease check the data & model specification")
}
}

plot.linkfunction  <- function(x,legend.loc="topright",...) UseMethod("plot.linkfunction")
