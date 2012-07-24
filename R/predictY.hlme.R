
predictY.hlme <- function(x,newdata,na.action=1,...){

if(missing(newdata)) stop("The argument newdata should be specified")
if(missing(x)) stop("The argument x should be specified")
if (!inherits(x, "hlme")) stop("use only with \"hlme\" objects")
# ad 2/04/2012 Xnames2
if (!all(x$Xnames2 %in% c(colnames(newdata),"intercept"))) {
cat("newdata should at least include the following covariates: ", "\n")
cat(x$Xnames2[-1], "\n")}
if (!all(x$Xnames2 %in% c(colnames(newdata),"intercept"))) stop("see above")
if (!inherits(newdata, "data.frame")) stop("newdata should be a data.frame object")


if(x$conv==1|x$conv==2) {



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

if(!(na.action%in%c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument") 

if(na.action==1){
	na.action=na.omit
}else{
	na.action=na.fail
}



### Traitement des donnees manquantes

mcall <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
#mcall <- x$call[c(1,match(c("data"),names(x$call),0))]
mcall$na.action <- na.action
mcall$data <- newdata1

# fixed
m <- mcall
m$formula <- formula(paste("~",x$call$fixed[3],sep=""))
m[[1]] <- as.name("model.frame")	
m <- eval(m, sys.parent()) 
na.fixed <- attr(m,"na.action")

# mixture
if(!is.null(x$call$mixture)){
	m <- mcall
	m$formula <- x$call$mixture
	m[[1]] <- as.name("model.frame")	
	m <- eval(m, sys.parent()) 
	na.mixture <- attr(m,"na.action")
}else{
	na.mixture <- NULL
}

# random
if(!is.null(x$call$random)){
	m <- mcall
	m$formula <- x$call$random
	m[[1]] <- as.name("model.frame")	
	m <- eval(m, sys.parent()) 
 	na.random <- attr(m,"na.action")
}else{
	na.random <- NULL
}
# classmb
if(!is.null(x$call$classmb)){ 
	m <- mcall	
	m$formula <- x$call$classmb	
	m[[1]] <- as.name("model.frame")	
	m <- eval(m, sys.parent()) 
 	na.classmb <- attr(m,"na.action")
}else{
	na.classmb <- NULL
}
## Table sans donnees manquante: newdata
na.action <- unique(c(na.fixed,na.mixture,na.random,na.classmb))
if(!is.null(na.action)){
	newdata1 <- newdata1[-na.action,]
}


## Construction de nouvelles var eplicatives sur la nouvelle table
## fixed
	
	X_fixed <- model.matrix(formula(paste("~",x$call$fixed[3],sep="")),data=newdata1)
	if(colnames(X_fixed)[1]=="(Intercept)"){
		colnames(X_fixed)[1] <- "intercept"
		int.fixed <- 1
	}	
## mixture
	if(!is.null(x$call$mixture)){
		X_mixture <- model.matrix(formula(x$call$mixture),data=newdata1)	
		if(colnames(X_mixture)[1]=="(Intercept)"){
			colnames(X_mixture)[1] <- "intercept"
			int.mixture <- 1
		}
		id.X_mixture <- 1
	}else{
		id.X_mixture <- 0
	}	
## random
	if(!is.null(x$call$random)){
		X_random <- model.matrix(formula(x$call$random),data=newdata1)	
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
		X_classmb <- model.matrix(formula(x$call$classmb),data=newdata1)
		colnames(X_classmb)[1] <- "intercept"
		id.X_classmb <- 1
	}else{
		id.X_classmb <- 0
	}
## Construction des var expli
newdata1 <- X_fixed


if(id.X_mixture == 1){
	for(i in 1:length(colnames(X_mixture))){
		if((colnames(X_mixture)[i] %in% colnames(newdata1))==F){
			newdata1 <- cbind(newdata1,X_mixture[,i])
			
		}
	}
}
if(id.X_random == 1){
	for(i in 1:length(colnames(X_random))){
		if((colnames(X_random)[i] %in% colnames(newdata1))==F){
			newdata1 <- cbind(newdata1,X_random[,i])
		}	 
	}
}
if(id.X_classmb == 1){
	for(i in 1:length(colnames(X_classmb))){
		if((colnames(X_classmb)[i] %in% colnames(newdata1))==F){
			newdata1 <- cbind(newdata1,X_classmb[,i],deparse.level=0)	 
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


if (x$ng==1){
colnames(Y) <- c("Ypred")
}
if (x$ng>1){
colnames(Y) <- c(paste("Ypred_class",1:x$ng,sep=""))
}
res <- Y
}else{
cat("Predictions can not be computed since the program stopped abnormally.")
}
}       


predictY <- function(x,...) UseMethod("predictY")
