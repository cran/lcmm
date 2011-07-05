

plot.predict.Jointlcmm <- function(x,newdata,var.time,legend.loc="topright",ylim=NULL,...){


if(missing(var.time)) stop("The argument var.time should be specified")
if(missing(newdata)) stop("The argument newdata should be specified")
if(missing(x)) stop("The argument x should be specified")
if (!inherits(x, "Jointlcmm")) stop("use only with \"Jointlcmm\" objects")
if(!all(x$Xnames %in% colnames(newdata))) {
cat("newdata should at least include the following covariates: ", "\n")
cat(x$Xnames, "\n")}
if (!all(x$Xnames %in% colnames(newdata))) stop("see above")
if (!inherits(newdata, "data.frame")) stop("newdata should be a data.frame object")
if (!inherits(var.time, "character")) stop("the class of var.time should be character")
if (!(var.time %in% colnames(newdata))) stop("The variable names var.time should be included in newdata")

if(x$conv==1|x$conv==2) {



newdata<-newdata[sort.list(newdata[,var.time]),] 
X <- newdata[,var.time]
newdata1 <- newdata[,x$Xnames]

X1 <- NULL
X2 <- NULL
b1 <- NULL
b2 <- NULL


kk<-0
for(k in 1:length(x$idg0)){
if(x$idg0[k]==1){
X1 <- cbind(X1,newdata1[,k])
place <- x$N[1]+x$N[2]+x$N[3]+kk
b1 <- c(b1,x$best[place+1])
kk <- kk+1}

if(x$idg0[k]==2){
X2 <- cbind(X2,newdata1[,k])
place1 <- x$N[1]+x$N[2]+x$N[3]+kk+1
place2 <- x$N[1]+x$N[2]+x$N[3]+kk+x$ng
b2 <- rbind(b2,x$best[place1:place2])
kk <- kk+x$ng}
}

Y<-matrix(0,length(newdata1[,1]),x$ng)
colnames(Y) <- paste("class",1:x$ng,sep="") 
for(g in 1:x$ng){
if(length(b1) != 0){
Y[,g]<- X1 %*% b1 
}
if(length(b2) != 0){
Y[,g]<- Y[,g] + X2 %*% b2[,g]
}
}

ylim1 <- ylim
if (is.null(ylim)){
ylim1 <- c(min(Y),max(Y))
}



title1 <- "Class-specific mean predicted trajectory"
for (i in 1:x$ng){
plot(Y[,i]~X,col=i,type="l",ylim=ylim1,xlim=c(min(X),max(X)),xlab=var.time,pch=2,ylab="outcome",main=title1,lty=i)
par(new=TRUE)
}
legend(x=legend.loc,legend=colnames(Y),col=1:x$ng,lty=1:x$ng,inset=.02)
par(new=FALSE)
}else{
cat("Predictions can not be computed since the program stopped abnormally.")
}
}


plot.predict <- function(x,...) UseMethod("plot.predict")
