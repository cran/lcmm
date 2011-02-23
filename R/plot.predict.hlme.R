

plot.predict.hlme <- function(x,newdata,var.time,legend.loc="topright",...){


if(missing(var.time)) stop("The argument var.time should be specified")
if(missing(newdata)) stop("The argument newdata should be specified")
if(missing(x)) stop("The argument x should be specified")
if (!inherits(x, "hlme")) stop("use only with \"hlme\" objects")
if(!identical(colnames(newdata),x$Xnames)) {
cat("newdata should include the following covariates in exactly the same order: ", "\n")
cat(x$Xnames, "\n")}
if (!identical(colnames(newdata),x$Xnames)) stop("see above")
if (!inherits(newdata, "data.frame")) stop("newdata should be a data.frame object")
if (!inherits(var.time, "character")) stop("the class of var.time should be character")



newdata<-newdata[sort.list(newdata[,var.time]),] 
X <- newdata[,var.time]


X1 <- NULL
X2 <- NULL
b1 <- NULL
b2 <- NULL

kk<-0
for(k in 1:length(x$idg0)){
if(x$idg0[k]==1){
X1 <- cbind(X1,newdata[,k])
place <- x$N[1]+kk
b1 <- c(b1,x$best[place+1])
kk <- kk+1}

if(x$idg0[k]==2){
X2 <- cbind(X2,newdata[,k])
place1 <- x$N[1]+kk+1
place2 <- x$N[1]+kk+x$ng
b2 <- rbind(b2,x$best[place1:place2])
kk <- kk+x$ng}
}

Y<-matrix(0,length(newdata[,1]),x$ng)
colnames(Y) <- paste("class",1:x$ng,sep="") 
for(g in 1:x$ng){
if(length(b1) != 0){
Y[,g]<- X1 %*% b1 
}
if(length(b2) != 0){
Y[,g]<- Y[,g] + X2 %*% b2[,g]
}
}

title1 <- "Class-specific mean predicted trajectory"
for (i in 1:x$ng){
plot(Y[,i]~X,col=i,type="l",ylim=c(min(Y),max(Y)),xlim=c(min(X),max(X)),xlab=var.time,pch=2,ylab="outcome",main=title1,lty=i)
par(new=TRUE)
}
legend(x=legend.loc,legend=colnames(Y),col=1:x$ng,lty=1:x$ng,inset=.02)
par(new=FALSE)
}


plot.predict <- function(x,...) UseMethod("plot.predict")
