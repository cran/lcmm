

predictY.Jointlcmm <- function(x,newdata,na.action=1,...){

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


if(x$conv==1|x$conv==2) {

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

call_fixed <- x$call$fixed[3]
if(is.null(x$call$random)) {call_random <- ~-1} else call_random <- x$call$random
if(is.null(x$call$classmb)) {call_classmb <- ~-1} else call_classmb <- x$call$classmb
if(is.null(x$call$mixture)) {call_mixture <- ~-1} else call_mixture <- x$call$mixture
if(is.null(x$call$survival)) {call_survival <- ~-1} else call_survival <- x$call$survival[3]


if(!(na.action%in%c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument") 

if(na.action==1){
	na.action=na.omit
}else{
	na.action=na.fail
}

### pour les facteurs

Xnames2 <- x$Names2[[1]]

 #cas où une variable du dataset est un facteur
 olddata <- eval(x$call$data)
  for(v in Xnames2[-1])
 {
  if (is.factor(olddata[,v]))     
  {
   mod <- levels(olddata[,v])
   if (!(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
   newdata1[,v] <- factor(newdata1[,v], levels=mod)
  }
 }
 
 #cas où on a factor() dans l'appel
 z <- all.names(call_fixed)
 ind_factor <- which(z=="factor")
 if(length(ind_factor))
 {
  nom.factor <- z[ind_factor+1]  
  for (v in nom.factor)
  {
   mod <- levels(as.factor(olddata[,v]))
   if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
   newdata1[,v] <- factor(newdata1[,v], levels=mod)
  }
 }
 call_fixed <- gsub("factor","",call_fixed)

 z <- all.names(call_random)
 ind_factor <- which(z=="factor")
 if(length(ind_factor))
 {
  nom.factor <- z[ind_factor+1]
  for (v in nom.factor)
  {
   mod <- levels(as.factor(olddata[,v]))
   if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
   newdata1[,v] <- factor(newdata1[,v], levels=mod)
  }
 }
 call_random <- gsub("factor","",call_random)
       
 z <- all.names(call_classmb)
 ind_factor <- which(z=="factor")
 if(length(ind_factor))
 {
  nom.factor <- z[ind_factor+1]
  for (v in nom.factor)
  {
   mod <- levels(as.factor(olddata[,v]))
   if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
   newdata1[,v] <- factor(newdata1[,v], levels=mod)
  }
 }
 call_classmb <- gsub("factor","",call_classmb)
        
 z <- all.names(call_mixture)
 ind_factor <- which(z=="factor")
 if(length(ind_factor))
 {
  nom.factor <- z[ind_factor+1]
  for (v in nom.factor)
  {
   mod <- levels(as.factor(olddata[,v]))
   if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
   newdata1[,v] <- factor(newdata1[,v], levels=mod)
  }
 }
 call_mixture <- gsub("factor","",call_mixture) 

 z <- all.names(call_survival)
 ind_factor <- which(z=="factor")
 if(length(ind_factor))
 {
  nom.factor <- z[ind_factor+1]  
  for (v in nom.factor)
  {
   mod <- levels(as.factor(olddata[,v]))
   if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
   newdata1[,v] <- factor(newdata1[,v], levels=mod)
  }
 }
 call_survival <- gsub("factor","",call_survival)
 call_survival <- gsub("mixture","",call_survival)
 
   
     
### Traitement des donnees manquantes
# fixed
#mcall <- x$call[c(1,match(c("data"),names(x$call),0))]
mcall <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
mcall$na.action <- na.action
mcall$data <- newdata1


# fixed
m <- mcall
m$formula <- formula(paste("~",call_fixed,sep=""))
m[[1]] <- as.name("model.frame")	
m <- eval(m, sys.parent()) 
na.fixed <- attr(m,"na.action")

# mixture
if(!is.null(x$call$mixture)){
	m <- mcall
	m$formula <- formula(paste("~",call_mixture,sep=""))
	m[[1]] <- as.name("model.frame")	
	m <- eval(m, sys.parent()) 
	na.mixture <- attr(m,"na.action")
}else{
	na.mixture <- NULL
}

# random
if(!is.null(x$call$random)){
	m <- mcall
	m$formula <- formula(paste("~",call_random,sep=""))
	m[[1]] <- as.name("model.frame")	
	m <- eval(m, sys.parent()) 
 	na.random <- attr(m,"na.action")
}else{
	na.random <- NULL
}

# classmb
if(!is.null(x$call$classmb)){ 
	m <- mcall	
	m$formula <- formula(paste("~",call_classmb,sep=""))
	m[[1]] <- as.name("model.frame")	
	m <- eval(m, sys.parent()) 
 	na.classmb <- attr(m,"na.action")
}else{
	na.classmb <- NULL
}

#survival
if(!is.null(x$call$survival))
{
 m <- mcall
 m$formula <- formula(paste("~",call_survival,sep=""))
 m[[1]] <- as.name("model.frame")	
 m <- eval(m, sys.parent()) 
 na.survival <- attr(m,"na.action")
}
else {na.survival <- NULL}

########### For survival
#if(!is.null(x$call$survival)){ 
#	res.evt <- x$call$survival
#	class(res.evt) <- "formula"
#	tmp.res <- res.evt
#
#	res.evt <- terms(res.evt,"mixture") 
#	inddep.surv <- attr(res.evt, "term.labels")
#	ind.mixture <- untangle.specials(res.evt, "mixture", 1)
#	inddepvar.Mixt  <- gsub("\\)","",gsub("mixture\\(","",ind.mixture$vars))
#	inddepvar.noMixt <- inddep.surv[!(inddep.surv %in% ind.mixture$vars)]
#	tmp.su <- c(inddepvar.noMixt,inddepvar.Mixt)
#	names.survival <- "~" 
#	for(i in 1:length(tmp.su)){
#	      if(i==1){
#		      names.survival <- paste(names.survival,tmp.su[i],sep="")
#	      }else{
#		      names.survival <- paste(names.survival,tmp.su[i],sep="+")
#	      }
#	}

##7/05/2012
##7/05/2012
#	m <- mcall
#	m$formula <- formula(names.survival)
#	m[[1]] <- as.name("model.frame")
#	m <- eval(m, sys.parent())
#	na.survival <- attr(m,"na.action")

## Table sans donnees manquante: newdata
	na.action <- unique(c(na.fixed,na.mixture,na.random,na.classmb,na.survival))
	if(!is.null(na.action)){
		newdata1 <- newdata1[-na.action,]
	}
##7/05/2012
#	X_survival <- model.matrix(formula(names.survival),data=newdata1)
#	if(colnames(X_survival)[1]=="(Intercept)"){
#		colnames(X_survival)[1] <- "intercept"
#	}
#	id.X_survival <- 1
#}else{
#	id.X_survival <- 0
##7/05/2012
#	na.action <- unique(c(na.fixed,na.mixture,na.random,na.classmb))
#	if(!is.null(na.action)){
#		newdata1 <- data[-na.action,]
#	}
#}



## Construction de nouvelles var explicatives sur la nouvelle table
## fixed
	
	X_fixed <- model.matrix(formula(paste("~",call_fixed,sep="")),data=newdata1)
	if(colnames(X_fixed)[1]=="(Intercept)"){
		colnames(X_fixed)[1] <- "intercept"
		int.fixed <- 1
	}	
## mixture
	if(!is.null(x$call$mixture)){
		X_mixture <- model.matrix(formula(paste("~",call_mixture,sep="")),data=newdata1)	
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
		X_random <- model.matrix(formula(paste("~",call_random,sep="")),data=newdata1)	
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
		X_classmb <- model.matrix(formula(paste("~",call_classmb,sep="")),data=newdata1)
		colnames(X_classmb)[1] <- "intercept"
		id.X_classmb <- 1
	}else{
		id.X_classmb <- 0
	}	

## survival
	if(!is.null(x$call$survival)){ 
		X_survival <- model.matrix(formula(paste("~",call_survival,sep="")),data=newdata1)
		colnames(X_survival)[1] <- "intercept"
		id.X_survival <- 1
	}else{
		id.X_survival <- 0
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


if (x$specif[[3]]==1){
colnames(Y) <- c("Ypred")
}
if (x$specif[[3]]>1){
colnames(Y) <- c(paste("Ypred_class",1:x$specif[[3]],sep=""))
}
res <- Y

}else{
cat("Predictions can not be computed since the program stopped abnormally.")
}
}


predictY <- function(x,...) UseMethod("predictY")
