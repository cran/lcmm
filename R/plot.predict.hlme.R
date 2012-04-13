

plot.predict.hlme <- function(x,newdata,var.time,legend.loc="topright",ylim=NULL,...){

if(missing(var.time)) stop("The argument var.time should be specified")
if(missing(newdata)) stop("The argument newdata should be specified")
if(missing(x)) stop("The argument x should be specified")
if (!inherits(x, "hlme")) stop("use only with \"hlme\" objects")
# ad 2/04/2012 Xnames2
if (!all(x$Xnames2 %in% c(colnames(newdata),"intercept"))) {
cat("newdata should at least include the following covariates: ", "\n")
cat(x$Xnames2[-1], "\n")}
if (!all(x$Xnames2 %in% c(colnames(newdata),"intercept"))) stop("see above")
if (!inherits(newdata, "data.frame")) stop("newdata should be a data.frame object")
if (!inherits(var.time, "character")) stop("the class of var.time should be character")
if (!(var.time %in% colnames(newdata))) stop("The variable names var.time should be included in newdata")

if(x$conv==1|x$conv==2) {

dev.new()


newdata<-newdata[sort.list(newdata[,var.time]),] 
X <- newdata[,var.time]

#------------> changement Cecile 10/04/2012
## add 12/04/2012
if(x$Xnames2[1]!="intercept"){
	newdata1 <- newdata[,x$Xnames2]
	colnames(newdata1) <- x$Xnames
	newdata1 <- data.frame(newdata1)
}else{
	newdata1 <- cbind(rep(1,length=length(newdata[,1])),newdata[,x$Xnames2[-1]])
	colnames(newdata1) <- c("intercept",x$Xnames2[-1])
	newdata1 <- data.frame(newdata1)
}



X1 <- NULL
X2 <- NULL
b1 <- NULL
b2 <- NULL


### add 11/04/2012
ff <- x$call
ff$fixed <- ff$mixture <- ff$random<- ff$subject<- ff$classmb<- ff$ng<- ff$idiag<-
ff$nwg<- ff$B<- ff$convB<- ff$convL<- ff$convG<- ff$prior<- ff$maxiter<-NULL


## fixed
fit <- ff
temp <- formula(x$call$fixed)
fit$formula <- formula(paste("~",temp[3],sep=""))
fit$data <- newdata1
fit[[1]] <- as.name("model.frame")
fit <- eval(fit, sys.parent()) 
X_fixed <- if (!is.empty.model(attr(fit, "terms")))model.matrix(attr(fit, "terms"), fit, contrasts)
if(colnames(X_fixed)[1]=="(Intercept)"){
	colnames(X_fixed)[1] <- "intercept"
}else{
	X_fixed <- cbind(intercept=rep(1,length(X_fixed[,1])),X_fixed)
}

## mixture
if(!is.null(x$call$mixture)){
	fit <- ff
	fit$formula <- x$call$mixture
	fit$data <- newdata1
	fit[[1]] <- as.name("model.frame")
	fit <- eval(fit, sys.parent()) 
	X_mixture <- if (!is.empty.model(attr(fit, "terms")))model.matrix(attr(fit, "terms"), fit, contrasts)
	if(colnames(X_mixture)[1]=="(Intercept)"){
		colnames(X_mixture)[1] <- "intercept"
		int.mixture <- 1
	}
	id.X_mixture <- 1
}else{
	id.X_mixture <- 0
}

## mixture
if(!is.null(x$call$random)){
	fit <- ff
	fit$formula <- x$call$random
	fit$data <- newdata1
	fit[[1]] <- as.name("model.frame")
	fit <- eval(fit, sys.parent()) 
	X_random <- if (!is.empty.model(attr(fit, "terms")))model.matrix(attr(fit, "terms"), fit, contrasts)
	if(colnames(X_random)[1]=="(Intercept)"){
		colnames(X_random)[1] <- "intercept"
		int.random <- 1
	}
	id.X_random <- 1
}else{
	id.X_random <- 0
}

## classmb
if(!is.null(x$call$classmb)){
	fit <- ff
	fit$formula <- x$call$classmb
	fit$data <- newdata1
	fit[[1]] <- as.name("model.frame")
	fit <- eval(fit, sys.parent()) 
	X_classmb <- if (!is.empty.model(attr(fit, "terms")))model.matrix(attr(fit, "terms"), fit, contrasts)
	if(colnames(X_classmb)[1]=="(Intercept)"){
		colnames(X_classmb)[1] <- "intercept"
		id.X_classmb <- 1
	}
	id.X_classmb <- 1
}else{
	id.X_classmb <- 0
}

newdata1 <- X_fixed
if(id.X_mixture == 1){
	for(i in 1:length(colnames(X_mixture))){
		if((colnames(X_mixture)[i] %in% colnames(newdata1))==F){
			newdata1 <- unique(cbind(newdata1,X_mixture[,i]))
			
		}
	}
}
if(id.X_random == 1){
	for(i in 1:length(colnames(X_random))){
		if((colnames(X_random)[i] %in% colnames(newdata1))==F){
			newdata1 <- unique(cbind(newdata1,X_random[,i]))
		}	 
	}
}
if(id.X_classmb == 1){
	for(i in 1:length(colnames(X_classmb))){
		if((colnames(X_classmb)[i] %in% colnames(newdata1))==F){
			newdata1 <- unique(cbind(newdata1,X_classmb[,i],deparse.level=0))	 
		}	
	}
}

##end add 11/04/2012

kk<-0
for(k in 1:length(x$idg0)){
if(x$idg0[k]==1){
X1 <- cbind(X1,newdata1[,k])
place <- x$N[1]+kk
b1 <- c(b1,x$best[place+1])
kk <- kk+1}

if(x$idg0[k]==2){
X2 <- cbind(X2,newdata1[,k])
place1 <- x$N[1]+kk+1
place2 <- x$N[1]+kk+x$ng
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

title1 <- "Class-specific mean latent process predicted trajectory"
for (i in 1:x$ng){
plot(Y[,i]~X,col=i,type="l",ylim=ylim1,xlim=c(min(X),max(X)),xlab=var.time,pch=2,ylab="outcome",main=title1,lty=i,...)
par(new=TRUE)
}
legend(x=legend.loc,legend=colnames(Y),col=1:x$ng,lty=1:x$ng,inset=.02,bty="n")
par(new=FALSE)
}else{
cat("Output can not be produced since the program stopped abnormally.")
}
}            


plot.predict <- function(x,...) UseMethod("plot.predict")
