

plot.predict.Jointlcmm <- function(x,newdata,var.time,legend.loc="topright",ylim=NULL,na.action=1,...){


if(missing(var.time)) stop("The argument var.time should be specified")
if(missing(newdata)) stop("The argument newdata should be specified")
if(missing(x)) stop("The argument x should be specified")
if (!inherits(x, "Jointlcmm")) stop("use only with \"Jointlcmm\" objects")
#--------------->
if(!all(x$Names2[[1]] %in% c(colnames(newdata),"intercept"))) {
cat("newdata should at least include the following covariates: ", "\n")
cat(x$Names2[[1]][-1], "\n")
}
if (!all(x$Names2[[1]] %in% c(colnames(newdata),"intercept"))) stop("see above")
#--------------->


if (!inherits(newdata, "data.frame")) stop("newdata should be a data.frame object")
if (!inherits(var.time, "character")) stop("the class of var.time should be character")
if (!(var.time %in% colnames(newdata))) stop("The variable names var.time should be included in newdata")

if(x$conv==1|x$conv==2) {

dev.new()


newdata<-newdata[sort.list(newdata[,var.time]),] 

#------------> changement Cecile 10/04/2012 ad
if(x$Names2[[1]][1]!="intercept"){
	newdata1 <- newdata[,x$Names2[[1]]]
	colnames(newdata1) <- x$Names2[[1]]
	newdata1 <- data.frame(newdata1)
}else{
	newdata1 <- cbind(intercept=rep(1,length=length(newdata[,1])),newdata[,x$Names2[[1]][-1]])
	colnames(newdata1) <- x$Names2[[1]]
	newdata1 <- data.frame(newdata1)
}
#------------>
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
# fixed
#mcall <- x$call[c(1,match(c("data"),names(x$call),0))]
mcall <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
mcall$na.action <- na.action
mcall$data <- newdata1


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

########## For survival
if(!is.null(x$call$survival)){ 
	res.evt <- x$call$survival
	class(res.evt) <- "formula"
	tmp.res <- res.evt

	res.evt <- terms(res.evt,"mixture") 
	inddep.surv <- attr(res.evt, "term.labels")
	ind.mixture <- untangle.specials(res.evt, "mixture", 1)
	inddepvar.Mixt  <- gsub("\\)","",gsub("mixture\\(","",ind.mixture$vars))
	inddepvar.noMixt <- inddep.surv[!(inddep.surv %in% ind.mixture$vars)]
	tmp.su <- c(inddepvar.noMixt,inddepvar.Mixt)
	names.survival <- "~" 
	for(i in 1:length(tmp.su)){
	      if(i==1){
		      names.survival <- paste(names.survival,tmp.su[i],sep="")
	      }else{
		      names.survival <- paste(names.survival,tmp.su[i],sep="+")
	      }
	}

#7/05/2012
#7/05/2012
	m <- mcall
	m$formula <- formula(names.survival)
	m[[1]] <- as.name("model.frame")
	m <- eval(m, sys.parent())
	na.survival <- attr(m,"na.action")

## Table sans donnees manquante: newdata
	na.action <- unique(c(na.fixed,na.mixture,na.random,na.classmb,na.survival))
	if(!is.null(na.action)){
		newdata1 <- newdata1[-na.action,]
	}
#7/05/2012
	X_survival <- model.matrix(formula(names.survival),data=newdata1)
	if(colnames(X_survival)[1]=="(Intercept)"){
		colnames(X_survival)[1] <- "intercept"
	}
	id.X_survival <- 1
}else{
	id.X_survival <- 0
#7/05/2012
	na.action <- unique(c(na.fixed,na.mixture,na.random,na.classmb))
	if(!is.null(na.action)){
		newdata1 <- data[-na.action,]
	}
}

X <- newdata1[,var.time]




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
if(id.X_survival == 1){
	for(i in 1:length(colnames(X_survival))){
		if((colnames(X_survival)[i] %in% colnames(newdata1))==F){
			newdata1 <- cbind(newdata1,X_survival[,i])
			
		}
	}
}

#colnames(newdata1) <- namesNew
### end add 11/04/2012



kk<-0
for(k in 1:length(x$specif[[6]])){
if(x$specif[[6]][k]==1){
X1 <- cbind(X1,newdata1[,k])
place <- x$specif[[1]][1]+x$specif[[1]][2]+x$specif[[1]][3]+kk
b1 <- c(b1,x$best[place+1])
kk <- kk+1}

if(x$specif[[6]][k]==2){
X2 <- cbind(X2,newdata1[,k])
place1 <- x$specif[[1]][1]+x$specif[[1]][2]+x$specif[[1]][3]+kk+1
place2 <- x$specif[[1]][1]+x$specif[[1]][2]+x$specif[[1]][3]+kk+x$specif[[3]]
b2 <- rbind(b2,x$best[place1:place2])
kk <- kk+x$specif[[3]]}
}

Y<-matrix(0,length(newdata1[,1]),x$specif[[3]])
colnames(Y) <- paste("class",1:x$specif[[3]],sep="") 
for(g in 1:x$specif[[3]]){
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
for (i in 1:x$specif[[3]]){
plot(Y[,i]~X,col=i,type="l",ylim=ylim1,xlim=c(min(X),max(X)),xlab=var.time,pch=2,ylab="outcome",main=title1,lty=i,...)
par(new=TRUE)
}
legend(x=legend.loc,legend=colnames(Y),col=1:x$specif[[3]],lty=1:x$specif[[3]],inset=.02,bty="n")
par(new=FALSE)
}else{
cat("Predictions can not be computed since the program stopped abnormally.")
}
}


plot.predict <- function(x,...) UseMethod("plot.predict")
