plot.linkfunction.lcmm <- function(x,legend.loc="topright",...){

if(missing(x)) stop("The argument x should be specified")
if (!inherits(x, "lcmm")) stop("use only with \"lcmm\" objects")
if(x$conv %in% c(1,2)){
dev.new()

if (x$linktype==3 & (x$linknodes[2]-x$linknodes[1])>1) {

ntrtot <- sum(x$ide==1)
diff <- x$estimlink[(2*(x$linknodes[2]-x$linknodes[1]+1)-1),2]-x$estimlink[2,2]
diff <- diff/ntrtot

lim <- c(x$estimlink[2,2]-diff,x$estimlink[(2*(x$linknodes[2]-x$linknodes[1]+1)-1),2]+diff)
}
else{
lim <- c(min(x$estimlink[,2]),max(x$estimlink[,2]))
}

title1 <- "Estimated link function"
plot(x$estimlink[,1]~x$estimlink[,2],type="l",xlim=lim,ylim=c(min(x$estimlink[,1]),max(x$estimlink[,1])),xlab="Latent process",pch=2,ylab="Longitudinal outcome",main=title1,...)
par(new=FALSE)
}else{
cat("Output can not be produced since the program stopped abnormally. \n")
}
}

plot.linkfunction  <- function(x,legend.loc="topright",...) UseMethod("plot.linkfunction")
