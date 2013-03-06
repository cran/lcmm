############## last change 2012/03/19 #####################

############################### Fonction Jointlcmm ###################################

Jointlcmm <-
function(fixed,mixture,random,subject,classmb,ng=1,idiag=FALSE,nwg=FALSE,survival,hazard="Weibull",hazardtype="Specific",hazardnodes=NULL,TimeDepVar=NULL,cor=NULL,data,B,convB=0.0001,convL=0.0001,convG=0.0001,maxiter=150,nsim=100,prior,logscale=FALSE,subset=NULL,na.action=1){

cl <- match.call()
m <- match.call(expand.dots = FALSE)
    m$fixed <- m$mixture <- m$random <- m$subject <- m$classmb <-m$ng<-m$idiag<-m$nwg<-m$B<-m$convB <- m$convL <-m$convG<-m$prior<-m$logscale<-m$maxiter<-m$hazardnodes<-m$hazard<-m$hazardtype<-m$TimeDepVar<-m$nsim<-m$... <-NULL
args <- as.list(match.call(Jointlcmm))[-1]

#subject.name <- as.character(args$subject)
subject.name <- as.character(subject)

#### ERROR MESSAGES
if(!missing(mixture) & ng==1) stop("No mixture can be specified with ng=1")
#if(missing(mixture) & ng>1) stop("The argument mixture is required for ng > 1")
if(!missing(classmb) & ng==1) stop("No classmb can be specified with ng=1")
if(missing(random)) random <- ~-1
if(missing(fixed)) stop("The argument Fixed must be specified in any model")
if(missing(classmb)) classmb <- ~-1
if(missing(mixture)) mixture <- ~-1
if(ng==1&nwg==TRUE) stop ("The argument nwg should be FALSE for ng=1")

if(class(fixed)!="formula") stop("The argument fixed must be a formula")
if(class(mixture)!="formula") stop("The argument mixture must be a formula")
if(class(random)!="formula") stop("The argument random must be a formula")
if(class(classmb)!="formula") stop("The argument classmb must be a formula")
if(class(survival)!="formula") stop("The argument survival must be a formula")
if(missing(data)){ stop("The argument data should be specified and defined as a data.frame")} 
if(nrow(data)==0) stop("Data should not be empty")
if(missing(subject)){ stop("The argument subject must be specified in any model even without random-effects")}

### test de l'argument cor
ncor0 <- 0
cor.type <- cl$cor[1]
cor.time <- cl$cor[2]
cor <- paste(cor.type,cor.time,sep="-")
if (all.equal(cor,character(0))!=TRUE)
{
 #if (typeof(cor)!= "character") stop("The argument cor must be a character")

 if (substr(cor,1,2)=="AR") { ncor0 <- 2 }
 else if (substr(cor,1,2)=="BM") { ncor0 <- 1  }
 else { stop("The argument cor must be of type AR or BM") }
 
 #if (substr(cor,3,3)!="-") stop("Invalid argument cor")
 #if (grep("-",cor)>1) stop("Invalid argument cor")
 
 if(!(strsplit(cor,"-")[[1]][2] %in% colnames(data))) stop("Unaible to find time variable from argument cor in data")
 else { cor.var.time <- strsplit(cor,"-")[[1]][2] }
}  
### fin test argument cor 

if(!(na.action%in%c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument") 

if(na.action==1){
	na.action=na.omit
}else{
	na.action=na.fail
}

### ad
X0.names2 <- c("intercept")
### ad
int.fixed <- 0
int.mixture <- 0
int.random <- 0
int.classmb <- 0

#7/05/2012
### Traitement des donnees manquantes
# fixed
m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
m$formula <- terms(fixed)
m$na.action <- na.action
m[[1]] <- as.name("model.frame")	
m <- eval(m, sys.parent()) 
na.fixed <- attr(m,"na.action")

# mixture
if(mixture[[2]] != "-1"){
	m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
	m$formula <- terms(mixture)
m$na.action <- na.action
	m[[1]] <- as.name("model.frame")	
	m <- eval(m, sys.parent()) 
	na.mixture <- attr(m,"na.action")	
}else{
	na.mixture <- NULL
}
# random
if(random[[2]] != "-1"){
	m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
	m$formula <- terms(random)
m$na.action <- na.action
	m[[1]] <- as.name("model.frame")	
	m <- eval(m, sys.parent()) 
 	na.random <- attr(m,"na.action")
}else{
	na.random <- NULL
}

# classmb
if(classmb[[2]] != "-1"){ 
	m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]	
	m$formula <- terms(classmb)
m$na.action <- na.action
	m[[1]] <- as.name("model.frame")	
	m <- eval(m, sys.parent()) 
 	na.classmb <- attr(m,"na.action")
}else{
	na.classmb <- NULL
}

#cor     
if(ncor0!=0)
{  
	m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
	m$formula <- as.formula(paste(cor.time,1,sep="~"))  
	m$na.action <- na.action
  m[[1]] <- as.name("model.frame")    
  m <- eval(m,sys.parent())  
  na.cor <- attr(m,"na.action") 	
}
else { na.cor <- NULL }


### names of covariate in intial fit
## fixed
X0.names2 <- unique(c(X0.names2,colnames(get_all_vars(formula(terms(fixed)),data=data))[-1]))
## mixture
if(mixture[[2]]!="-1")X0.names2 <- unique(c(X0.names2,colnames(get_all_vars(formula(terms(mixture)),data=data))))	
## random
if(random[[2]] != "-1"){
	X0.names2 <- unique(c(X0.names2,colnames(get_all_vars(formula(terms(random)),data=data))))
	inddepvar.random.nom2 <- colnames(get_all_vars(formula(terms(random)),data=data))
}
## classmb
if(classmb[[2]] != "-1")X0.names2 <- unique(c(X0.names2,colnames(get_all_vars(formula(terms(classmb)),data=data))))
#7/05/2012 

########################################################
#####              SURVIVAL  
########################################################
special <- c("mixture")

res.evt <- if (missing(data)){ 
		terms(survival, special)
	}else{
		terms(survival, special, data = data)  
	}  
	ord <- attr(res.evt, "order")
	if (length(ord) & any(ord != 1)){ 
		stop("Interaction terms are not valid for this function")
	}

### for covariates in survival
### nom des variables a droite dans survie (avec et sans mixture)
	inddep.surv <- attr(res.evt, "term.labels")

### ind.mixture = toutes les var avec mixture dans survie
	ind.mixture <- untangle.specials(res.evt, "mixture", 1)

### stop if mixture and ng=1 ad 10/04/2012

	if((ng == 1) && length(ind.mixture$terms)) stop("A class-specific effect of covariate can be specified in the survival model with mixture() only when number of groups is greater than 1. ")

### variables avec mixture
	dcomp <- ind.mixture$vars
	dcomp <- gsub("mixture\\(","",dcomp)
	dcomp <- gsub("\\)","",dcomp)
	inddepvar.Mixt  <- dcomp
## ad
	if(length(grep("factor",inddepvar.Mixt))!= 0)inddepvar.Mixt <- paste(inddepvar.Mixt,")",sep="")
## ad
### inddepvar.survival = toutes les variables avec pas de mixture dans survie
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
	int.survival <- 0
#7/05/2012
	m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
	m$formula <- formula(update(survival,names.survival))
m$na.action <- na.action
	m[[1]] <- as.name("model.frame")
	m <- eval(m, sys.parent())
	na.survival <- attr(m,"na.action")
#7/05/2012
	mtemp <- get_all_vars(formula(names.survival),data=data)
	X0.names2 <- unique(c(X0.names2,colnames(mtemp)))

#7/05/2012
## Table sans donnees manquante: newdata
	na.action <- unique(c(na.fixed,na.mixture,na.random,na.classmb,na.cor))
    # dans na.action, on a les indices des NA dans le subset de data

newdata <- data  
#prendre le subset :
if(!(is.null(subset))) 
newdata <- data[subset,]

#enlever les NA
	if(!is.null(na.action)){
		newdata <- newdata[-na.action,]
	}
#7/05/2012
## Construction de nouvelles var eplicatives sur la nouvelle table
## fixed
	X_fixed <- model.matrix(fixed,data=newdata)
	if(colnames(X_fixed)[1]=="(Intercept)"){
		colnames(X_fixed)[1] <- "intercept"
		int.fixed <- 1
	}
	nom.fixed <- inddepvar.fixed <- inddepvar.fixed.nom <- colnames(X_fixed)
	if(int.fixed>0)inddepvar.fixed <- inddepvar.fixed[-1]

## mixture
	if(mixture[[2]] != "-1"){
		X_mixture <- model.matrix(mixture,data=newdata)	
		if(colnames(X_mixture)[1]=="(Intercept)"){
			colnames(X_mixture)[1] <- "intercept"
			int.mixture <- 1
		}
		nom.mixture <- inddepvar.mixture <- inddepvar.mixture.nom <- colnames(X_mixture)
		if(int.mixture>0)inddepvar.mixture <- inddepvar.mixture[-1]
		id.X_mixture <- 1
	}else{
		inddepvar.mixture <- nom.mixture <- inddepvar.mixture.nom <- NULL
		id.X_mixture <- 0
	}
## random
	if(random[[2]] != "-1"){
		X_random <- model.matrix(random,data=newdata)	
		if(colnames(X_random)[1]=="(Intercept)"){
			colnames(X_random)[1] <- "intercept"
			int.random <- 1
		}
		inddepvar.random <- inddepvar.random.nom <- colnames(X_random)
		if(int.random>0) inddepvar.random <- inddepvar.random[-1]
		id.X_random <- 1
	}else{
	## ad: add inddepvar.random.nom2 <- NULL 10/04/2012
		inddepvar.random <- inddepvar.random.nom <- inddepvar.random.nom2 <- NULL
		id.X_random <- 0
	}
## classmb
	if(classmb[[2]] != "-1"){ 
		X_classmb <- model.matrix(classmb,data=newdata)
		colnames(X_classmb)[1] <- "intercept"
		id.X_classmb <- 1
		inddepvar.classmb <- colnames(X_classmb)[-1]
		inddepvar.classmb.nom <- colnames(X_classmb)
	}else{
		inddepvar.classmb <- inddepvar.classmb.nom <- "intercept"
		id.X_classmb <- 0
	}	


	X_survival <- model.matrix(formula(update(survival,names.survival)),data=newdata)
	if(colnames(X_survival)[1]=="(Intercept)"){
		colnames(X_survival)[1] <- "intercept"
		int.survival <- 1
	}

	if(length(inddepvar.Mixt) > 0){
		inddepvar.Mixt <- colnames(X_survival)[grep(inddepvar.Mixt,colnames(X_survival))]
		inddepvar.noMixt <-colnames(X_survival)[!(colnames(X_survival) %in% inddepvar.Mixt)][-1]
	}else{
		inddepvar.noMixt <- colnames(X_survival)[-1]
	}

### ind.mixture$vars et ajouter a inddepvar.survival
    if(int.survival>0) inddepvar.survival <- colnames(X_survival)[-1] 


################# type du modele de survie
	depvar1 <- model.response(model.frame(formula=update(survival,"~1"),data=newdata))

	if (!inherits(depvar1, "Surv")){ 
		stop("Response must be a survival object")
	}

	attr.surv.type <- attr(depvar1,"type")
 
	if (attr.surv.type=="right"){
		idtrunc0 <-0
		Tevent <- depvar1[,1]
		Devent <- depvar1[,2]
		Tentry <- rep(0,length(Tevent))
		X0.names2 <- unique(c(X0.names2,(colnames(mtemp)[-c(1,2)]))) ### ad
	}

    if (attr.surv.type=="counting"){
        idtrunc0 <- 1
        Tentry <- depvar1[,1]
        Tevent <- depvar1[,2]
        Devent <- depvar1[,3]
	X0.names2 <- unique(c(X0.names2,(colnames(mtemp)[-c(1,2,3)]))) ### ad
    }

    if (max(Devent)==1) evt <- 1
    if (max(Devent)==0) evt <- 0

if (!(attr.surv.type %in% c("right","counting"))){
   stop("Jointlcmm handles only 'right' and 'counting' types of Survival data")
}


##############   COVARIATES       ##########################

##### debut ajout AMADOU

var.exp <- NULL
var.exp <- c(var.exp,colnames(X_fixed))
if(id.X_mixture == 1) var.exp <- c(var.exp,colnames(X_mixture))
if(id.X_random == 1)var.exp <- c(var.exp,colnames(X_random))
if(id.X_classmb == 1)var.exp <- c(var.exp,colnames(X_classmb))
var.exp <- unique(var.exp)
if(ncor0!=0)
{ if(!(cor.var.time %in% var.exp)) 
  {var.exp <- c(var.exp, cor.var.time)} #si la varaible de temps dans cor n'est dan sles variables expl, on l'ajoute
}

# intercept is always in inddepvar.classmb
var.exp <- unique(c(var.exp,colnames(X_survival)))

if(!(all(nom.mixture %in% nom.fixed))) stop("The covariates included in mixture should also be included in the argument fixed")


  
##############   DATA      ##########################
Y.name <- as.character(attributes(terms(fixed))$variables[2])
Y0 <- newdata[,Y.name]


##### debut ajout AMADOU
X0 <- X_fixed

oldnames <- colnames(X0)
if(id.X_mixture == 1){
	for(i in 1:length(colnames(X_mixture))){
		if((colnames(X_mixture)[i] %in% colnames(X0))==F){
			X0 <- cbind(X0,X_mixture[,i])
			colnames(X0) <- c(oldnames, colnames(X_mixture)[i])
			oldnames <- colnames(X0)
		}
	}
}   
if(id.X_random == 1){
	for(i in 1:length(colnames(X_random))){
		if((colnames(X_random)[i] %in% colnames(X0))==F){
			X0 <- cbind(X0,X_random[,i])
		  colnames(X0) <- c(oldnames, colnames(X_random)[i])
			oldnames <- colnames(X0)
		}	 
	}
}      
  
if(id.X_classmb == 1){
	for(i in 1:length(colnames(X_classmb))){
		if((colnames(X_classmb)[i] %in% colnames(X0))==F){
			X0 <- cbind(X0,X_classmb[,i], deparse.level=1)
      colnames(X0) <- c(oldnames,colnames(X_classmb)[i])
      oldnames <- colnames(X0)	
		}	
	}
}    
   
for(i in 1:length(colnames(X_survival))){
	if((colnames(X_survival)[i] %in% colnames(X0))==F){
		X0 <- cbind(X0,X_survival[,i], deparse.level=1) 
		colnames(X0) <- c(oldnames, colnames(X_survival)[i])
		oldnames <- colnames(X0)
	}
}     
if(ncor0!=0)
{ if(!(cor.var.time %in% colnames(X0))) 
  {
   X0 <- cbind(X0, newdata[,cor.var.time])
   colnames(X0) <- c(oldnames, cor.var.time)
  }
}  

if (sum(colnames(X0)==var.exp)!=length(var.exp)) print("Erreur sur les variables explicatives")
#colnames(X0) <- var.exp
  
#X0 <- X0[,-which(colnames(X0)=="intercept")]

if((any(is.na(X0))==TRUE)|(any(is.na(Y0))==TRUE))stop("Jointlcmm has not correctly handled missing data. Please revise 'na.action' option or remove manually missing data in the dataset")
 
n <- dim(data)[1]
#if ((int.fixed+int.random)>0) X0<- cbind(intercept=rep(1,n),X0)
# ad 10/04/2012
if (!((int.fixed+int.random)>0)) X0 <- as.data.frame(X0[,-which(colnames(X0)=="intercept")])

X0.names <- colnames(X0)
nvar.exp <- length(X0.names)

IND <- newdata[,subject.name]
IDnum <- as.numeric(IND)


#### INCLUSION PRIOR 
if(missing(prior)){ 
PRIOR <- seq(0,length=length(IND))
prior.name <- NULL
} 
if(!missing(prior)){ 
prior.name <- as.character(args$prior)
PRIOR <- newdata[,prior.name]
PRIOR[(is.na(PRIOR))] <- 0
}


#### Parametre timeDepVar  ####


TimeDepVar.name <- as.character(TimeDepVar)
nvdepsurv  <- 0
Tint <- Tevent
if(!missing(TimeDepVar)) {
### test : si pas une variable de data : probleme
Tint <- newdata[,TimeDepVar.name]
Tint[(is.na(Tint))] <- Tevent[(is.na(Tint))]
Tint[Tint>Tevent] <- Tevent[Tint>Tevent]
Tint[Tint<Tentry] <- Tentry[Tint<Tentry]
nvdepsurv <- 1
if (length(Tint[Tint>Tentry])==0) {
Tint <- Tevent
cat("TimeDepVar will be ignored since it is always lower than Time of Entry (0 by default). \n")
nvdepsurv  <- 0
}
if (length(Tint[Tint<Tevent])==0) {
cat("TimeDepVar will be ignored since it is always greater than Time of Event. \n")
nvdepsurv <- 0
}
}





#### DEFINITION INDICATORS 
ng0 <- ng
idiag0 <- as.integer(idiag)
nwg0 <- as.integer(nwg)
logspecif <- as.integer(logscale)
#cat("specif",logspecif)
idea0 <- rep(0,nvar.exp)
idprob0 <- rep(0,nvar.exp)
idg0 <- rep(0,nvar.exp)
idcor0 <- rep(0,nvar.exp)
idxevt <- rep(0,nvar.exp)

for (i in 1:nvar.exp){
 idea0[i] <- X0.names[i]%in%inddepvar.random.nom
 idprob0[i] <- X0.names[i]%in%inddepvar.classmb.nom      
 if(X0.names[i]%in%nom.fixed & !(X0.names[i]%in%nom.mixture)) idg0[i] <- 1 
 if(X0.names[i]%in%nom.fixed & X0.names[i]%in%nom.mixture) idg0[i] <- 2 
 if((X0.names[i]%in%inddepvar.survival) & !(X0.names[i]%in%inddepvar.Mixt)) idxevt[i] <- 1
 if((X0.names[i]%in%inddepvar.survival) & (X0.names[i]%in%inddepvar.Mixt)) idxevt[i] <- 2 
 }

if (ncor0!=0) idcor0 <- as.numeric(X0.names %in% cor.var.time) 
 
if((int.fixed+int.random)>0) idprob0[1] <- 0


###### DATA SORTING on IND variables
matYX <- cbind(IDnum,IND,PRIOR,Tentry,Tevent,Devent,Tint,Y0,X0)
matYXord <- matYX[sort.list(matYX[,1]),]

Y0 <- matYXord[,8]
X0 <- matYXord[,-c(1,2,3,4,5,6,7,8)]
IDnum <- matYXord[,1]
IND <-  matYXord[,2]
PRIOR <- matYXord[,3]
PRIOR <-as.integer(as.vector(PRIOR))
Tevent <- matYXord[,5]
Tentry <- matYXord[,4]
Devent <- matYXord[,6]
Tint <- matYXord[,7]
Devent<-as.integer(as.matrix(Devent))
Tevent<-as.numeric(as.matrix(Tevent))
Tentry<-as.numeric(as.matrix(Tentry))
Tint<-as.numeric(as.matrix(Tint))

X0<-as.numeric(as.matrix(X0))
Y0<-as.numeric(as.matrix(Y0))
nmes0<-as.vector(table(IDnum))
ns0<-length(nmes0)


##### INCLUSION PRIOR 
# definition de prior a 0 pour l'analyse G=1
initial <- as.integer(rep(0,ns0))
prior0 <- initial
# si prior pas missing alors mettre dedans la classe a priori. Attention tester q les valeurs sont dans 0, G
if(!missing(prior)){ 
prior0 <- PRIOR[cumsum(nmes0)]
}
INDuniq <- IND[cumsum(nmes0)]
seqnG <- 0:ng0
if (!(all(prior0  %in% seqnG))) stop ("The argument prior should contain integers between 0 and ng")


### reduction de Devent, Tevent,Tentry a la taille ns0
devt <- initial
tsurv0 <- initial
tsurv <- initial
devt <- Devent[cumsum(nmes0)]
tsurv <- Tevent[cumsum(nmes0)]
tsurv0 <- Tentry[cumsum(nmes0)]
tsurvint <- initial
tsurvint <- Tint[cumsum(nmes0)]
indsurvint <- initial
indsurvint[tsurvint<tsurv] <- 1


### variable  initialisation 

loglik <- as.double(0)
ni <- 0
istop <- 0
statsc <- 0
gconv <-rep(0,3)
ppi0 <- rep(0,ns0*ng0)
ppitest0 <- rep(0,ns0*ng0)
nv0<-nvar.exp
nobs0<-length(Y0)
resid_m <- rep(0,nobs0)
resid_ss <- rep(0,nobs0)
pred_m_g <- rep(0,nobs0*ng0)
pred_ss_g <- rep(0,nobs0*ng0)
nea0 <- sum(idea0==1)
predRE <- rep(0,nea0*ns0)
risq_est <- rep(0,nsim*ng0)
time_est <- rep(0,nsim)
surv_est <- rep(0,nsim*ng0)

nrisqtot <- 0
zitype <- 0
# zitype =1 for quant, 2 for equi, 3 for manual



##### hazard specification #########################################

# revoir le test ci-dessous : peut-�tre inclu dans un autre
carac <- strsplit(hazard,split="-")
carac <- unlist(carac)   
if(!(length(carac) %in% c(1,2,3))){stop("Please check and revise the hazard argument according to the format specified in the help.")}
     
# test sur hazardtype
if (!(hazardtype %in% c("PH","Common","Specific"))) stop("Only 'Specific', 'PH' or 'Common' hazardtype can be specified corresponding respectively to class-specific hazards, hazards proportional over latent classes or hazards common over classes")
     
haz <-strsplit(hazard,split="")
haz <- unlist(haz)
haz <- grep("-",haz)

### longueur hazard = 1
if((all.equal(length(haz),0)==T)==T){
   if(!(hazard %in% c("Weibull","piecewise","splines"))){
	stop("Only 'Weibull', 'piecewise' or 'splines' hazard can be specified in hazard argument")
   }
   else {
	nz0 <- switch(hazard,"Weibull"=2,"piecewise"=5,"splines"=5)
	typrisq0 <- switch(hazard,"Weibull"=2,"piecewise"=1,"splines"=3)
	nprisq0 <- switch(hazard,"Weibull"=2,"piecewise"=nz0-1,"splines"=nz0+2)
	zitype<- switch(hazard,"Weibull"=0,"piecewise"=1,"splines"=1)
   } 
}else{
#### longueur hazard > 1

 ## controle longueur
   if(!(length(haz) %in% c(1,2))) stop("With splines or piecewise baseline function, the separator of hazard argument must be only '-'")  
	
   test <- strsplit(hazard,split="-")
   test <- unlist(test)   
	
   if(all.equal(length(test),2)==T){
   ## si hazard a 2 elements
       if(!all(test[1:2] %in% c("splines","piecewise","equi","manual","quant"))){
		stop ("The hazard argument is incorrectly specified. Please refer to the help file of Joinlcmm.")
	 }
 	 # manuel
	 if(("manual" %in% unlist(strsplit(hazard,split="-")))){	
	      zitype <- 3 
		nz0 <- length(hazardnodes)+2
		typrisq0 <- switch(test[2],"piecewise"=1,"splines"=3)
		nprisq0 <- switch(test[2],"piecewise"=nz0-1,"splines"=nz0+2)
	  }
        else {
      	 # quant      
      	 if(("quant" %in% unlist(strsplit(hazard,split="-")))){	zitype <- 1 }
		 # equi
 		 if(("equi" %in% unlist(strsplit(hazard,split="-")))){	zitype <- 2 }
	       nz0 <- 5
		 typrisq0 <- switch(test[2],"piecewise"=1,"splines"=3)
		 nprisq0 <- switch(test[2],"piecewise"=nz0-1,"splines"=nz0+2)
	 } 
		
   }else{
   ## si hazard a 3 elements
       if(!all(test[2:3] %in% c("splines","piecewise","equi","manual","quant"))){
		stop ("The hazard argument is incorrectly specified. Please refer to the help file of Joinlcmm.")
	  }
        # manuel
	  if(("manual" %in% unlist(strsplit(hazard,split="-")))){	
	      zitype <- 3 
            nz0 <- as.integer(test[1])
		nz1 <- length(hazardnodes)+2
		if(!(all.equal(nz0,nz1)==T)){
			cat("Warning: the number of internal nodes does not correspond to the number of nodes indicated in 'hazard'.","\n")
			cat("The number of nodes derived from the list of internal nodes is kept for the analysis","\n") 
			nz0 <- nz1				
		}		
		typrisq0 <- switch(test[3],"piecewise"=1,"splines"=3)
		nprisq0 <- switch(test[3],"piecewise"=nz0-1,"splines"=nz0+2)
	  }
        else {
	       # quant      
      	 if(("quant" %in% unlist(strsplit(hazard,split="-")))){	zitype <- 1 }
		 # equi
 		 if(("equi" %in% unlist(strsplit(hazard,split="-")))){	zitype <- 2 }
	  	 nz0 <- as.integer(test[1])
		 typrisq0 <- switch(test[3],"piecewise"=1,"splines"=3)
		 nprisq0 <- switch(test[3],"piecewise"=nz0-1,"splines"=nz0+2)
	  }
   }	 
  
}

risqcom0 <- switch(hazardtype,"Specific"=0,"PH"=2,"Common"=1)
if (evt != 0){   
  nrisqtot <-  switch(hazardtype,"Specific"=nprisq0*ng0,"PH"=nprisq0+ng0-1,"Common"=nprisq0)   
} 


##### parametre nvarxevt  	    
nvarxevt <- sum(idxevt==1)+(sum(idxevt==2))*ng0 + nvdepsurv
if (evt==0) nvarxevt <- 0


##########   zi0 #################
#### Calcul de zi sur tsurv pour devt=1 ####
zi0 <- rep(0,nz0)
#TSURV tsurv pour devt=1
ns1 <- sum(devt==1)
TSURV <- rep(0,ns1)
k <- 0
for(i in 1:ns0){
   if(devt[i]==1){
      k <- k+1
      TSURV[k] <- tsurv[i]
   }
}
if(idtrunc0==1){
zi0[1] <- min(tsurv0)
zi0[nz0] <- max(c(tsurv0,tsurv))
}
if(idtrunc0==0){
#zi0[1] <- min(c(tsurv),c(tsurvint))
zi0[1] <- 0
zi0[nz0] <- max(tsurv)
}
minT <- zi0[1]
maxT <- zi0[nz0]
if((maxT-minT) < 0) stop("Please check the time of event variable. It seems that all the times are equal.")

##### si noeuds a initialiser
if (zitype > 0){

    if(nz0-2 <= 0){
         stop("Splines or piecewise baseline function should include at least 2 nodes (and at least 5 are recommended)")
    }
# equi
    if(zitype == 2){
         pas=as.double(maxT-minT)/as.double(nz0-1)
         for(i in 2:(nz0-1)){
              zi0[i] <- zi0[i-1]+ pas
         }
     }
#  manual
     if(zitype == 3){
        if (is.null(hazardnodes)){
             stop("If 'manual' option is specified for the splines or piecewise baseline hazard function, hazardnodes argument should include the list of interior nodes")
        }else{
             hazardnodes <- sort(hazardnodes)
             zi0[2:(nz0-1)] <- hazardnodes[1:(nz0-2)]
         }
     }
#  quant
     if(zitype == 1){
        pas <-c(1:(nz0-2))/(nz0-1) 
        zi0 [2:(nz0-1)] <- quantile(TSURV,probs=pas)
     }
}







#-------------------------------------------------------------------------------
#definition du vecteur de parametre + initialisation
#-------------------------------------------------------------------------------
#####cas 1 : ng=1


b<-NULL
b1 <- NULL
NPROB <- 0
###############################################################################

if(ng0==1| missing(B)){
   NVARXEVTinit <- sum(idxevt!=0) + nvdepsurv

   if (typrisq0==2) {
      b1[1] <- log(sum(devt)/sum(tsurv[devt==1]))
      b1[2] <- 0
   } 
   else {
     b1[1:(nprisq0)]<-(1/nprisq0)
   }
   b1[(1+nprisq0):(nprisq0+NVARXEVTinit)]<-0
   NEF<-nprisq0+NVARXEVTinit+sum(idg0!=0)
   b1[(nprisq0+NVARXEVTinit+1):(NEF)]<-0
   if(int.fixed > 0)  b1[(nprisq0+NVARXEVTinit+1)]<-mean(Y0)


   if(idiag0==1){
      NVC<-sum(idea0==1)
      b1[(NEF+1):(NEF+NVC)]<-1
   }
      
   if(idiag0==0){
      kk<-sum(idea0==1) 
      NVC<-(kk*(kk+1))/2
      indice<-cumsum(1:kk)
      bidiag<-rep(0,NVC)
      bidiag[indice]<-1
      b1[(NEF+1):(NEF+NVC)]<-bidiag
   }
   if(ncor0==1)
   {b1[NEF+NVC+1] <- 1 }
   if(ncor0==2)
   {b1[(NEF+NVC+1):(NEF+NVC+ncor0)] <- c(1,0) } 
   b1[(NEF+NVC+ncor0+1)]<-1
   NPM<-length(b1)
   NW<-0
   V <- rep(0,NPM*(NPM+1)/2) 
}


###############################################################################
#####cas 2 : ng>=2
if(ng0>1){
      NPROB<-(sum(idprob0==1)+1)*(ng0-1)
      b[1:NPROB]<-0
#      cat("nprob",NPROB,"nrisqtot",nrisqtot,"nvarxevt",nvarxevt,"\n")
      b[(NPROB+1):(NPROB+nrisqtot)]<- 1/nprisq0   
      b[(NPROB+nrisqtot+1):(NPROB+nrisqtot+nvarxevt)]<-0   
      
      NEF<-NPROB+nrisqtot+nvarxevt+sum(idg0==1)+(sum(idg0==2))*ng0      
      b[(NPROB+nrisqtot+nvarxevt+1):NEF]<-0
#      cat("nef",NEF,"\n")
      if(idiag0==1)NVC<-sum(idea0==1)
      if(idiag0==0){
         kk<-sum(idea0==1) 
         NVC<-(kk*(kk+1))/2
      }
#      cat("NVC",NVC,"\n")
      b[(NEF+1):(NEF+NVC)]<-1
      NW<-nwg0*(ng0-1)
#      cat("NW",NW,"\n")
      if(NW>0) b[(NEF+NVC+1):(NEF+NVC+NW)]<-1
      if(ncor0==1)
      {b[NEF+NVC+NW+1] <- 1 }
      if(ncor0==2)
      {b[(NPROB+NEF+NVC+NW+1):(NPROB+NEF+NVC+NW+ncor0)] <- c(0,1) }
      b[(NEF+NVC+NW+ncor0+1)]<-1
      NPM<-NEF+NVC+NW+ncor0+1
      V <- rep(0,NPM*(NPM+1)/2)
} 

###############################################################################
if(missing(B)){
    if(ng0>1){
         idea2 <- idea0
         idprob2 <- rep(0,nv0)  
         idg2 <- rep(0,nv0) 
         idg2[idg0!=0] <- 1
         idxevt2 <- rep(0,nv0) 
         idxevt2[idxevt!=0] <- 1 
         NEF2<-nprisq0+NVARXEVTinit+sum(idg2==1)
         NPM2<-NEF2+NVC+ncor0+1
         nwg2<-0
         ng2<-1
         ppi2 <- rep(0,ns0)
	       ppitest2 <- rep(0,ns0)
         pred_m_g2 <- rep(0,nobs0)
         pred_ss_g2 <- rep(0,nobs0)
         V2 <- rep(0,NPM2*(NPM2+1)/2)
         maxiter2 <- min(75,maxiter)
         convB2 <- max(0.01,convB)
         convL2 <- max(0.01,convL)
         convG2 <- max(0.01,convG)

init <- .Fortran("Jointhet",as.double(Y0),as.double(X0),as.integer(initial),as.integer(idprob2),
as.integer(idea2),as.integer(idg2),as.integer(idcor0),as.integer(idxevt2),as.integer(ns0),
as.integer(ng2),as.integer(nv0),as.integer(nobs0),as.integer(nea0),as.integer(ncor0),
as.integer(nmes0),as.integer(idiag0),as.integer(nwg2),as.integer(NPM2),best=as.double(b1),
V=as.double(V2),as.double(loglik),as.integer(ni),conv=as.integer(istop),as.double(gconv),as.double(ppi2),
as.double(ppitest2),as.double(resid_m),as.double(resid_ss),as.double(pred_m_g2),
as.double(pred_ss_g2),as.double(predRE),as.double(convB2),as.double(convL2),as.double(convG2),
as.integer(maxiter2),as.integer(typrisq0),as.integer(idtrunc0),as.integer(risqcom0),as.integer(nz0),
as.double(zi0),as.double(tsurv0),as.double(tsurv),as.double(tsurvint),as.integer(devt),
as.integer(indsurvint),as.double(statsc),as.double(risq_est),as.double(surv_est),as.integer(nsim),
as.double(time_est),as.integer(logspecif),PACKAGE="lcmm")


#      cat("Best",init$best,"\n")


	if (risqcom0==1) b[(NPROB+1):(NPROB+nprisq0)] <- init$best[1:nprisq0]   
	if (risqcom0==2) {
b[(NPROB+1):(NPROB+nprisq0)] <- init$best[1:nprisq0]   
b[(NPROB+nprisq0+1):(NPROB+nprisq0+ng-1)] <- 1+ (1:(ng-1))/(ng-1) 
}
	if (risqcom0==0) {
              for (g in 1:ng){
	           ident<- 1:nprisq0	
             ident <- ident*(ident+1)/2	         
	           if(init$conv==1)b[(NPROB+nprisq0*(g-1)+1):(NPROB+nprisq0*g)] <- abs(init$best[1:nprisq0])+(g-(ng+1)/2)*sqrt(init$V[ident])
	           else b[(NPROB+nprisq0*(g-1)+1):(NPROB+nprisq0*g)] <- abs(init$best[1:nprisq0])+(g-(ng+1)/2)*init$best[1:nprisq0]
              }
}
        k <- NPROB+nrisqtot
        l <- nprisq0
        t<- 0
        for (i in 1:nvar.exp){
           if(idxevt[i]==1){
              l <- l+1
              t <- t+1
              b[k+t] <- init$best[l]
           }
           if(idxevt[i]==2){
              l <- l+1
              for (g in 1:ng){
	         t <- t+1
	         if(init$conv==1) b[k+t] <- init$best[l]+(g-(ng+1)/2)*sqrt(init$V[l*(l+1)/2])
	         else b[k+t] <- init$best[l]+(g-(ng+1)/2)*init$best[l]
              }
           }
        }

## k pour le nouveau vecteur B (ie nrisqtot+nvarxevt)
## l pour le vecteur B issu de init (ie nprisq0+NVARXEVTinit)

        k <- NPROB+nrisqtot+nvarxevt
	 if (nvdepsurv==1) b[k] <- 0
        l <- nprisq0+NVARXEVTinit
        t<- 0
        for (i in 1:nvar.exp){
           if(idg0[i]==1){
              l <- l+1
              t <- t+1
              b[k+t] <- init$best[l]
           }
           if(idg0[i]==2){
              l <- l+1
              for (g in 1:ng){
	         t <- t+1
	         if(init$conv==1) b[k+t] <- init$best[l]+(g-(ng+1)/2)*sqrt(init$V[l*(l+1)/2])
	         else b[k+t] <- init$best[l]+(g-(ng+1)/2)*init$best[l]
              }
           }
        }
        b[(NEF+1):(NEF+NVC)] <-init$best[(NEF2+1):(NEF2+NVC)]
        b[(NEF+NVC+NW+ncor0+1)] <- init$best[NPM2]
        if (ncor0>0) {b[(NEF+NVC+NW+1):(NEF+NVC+NW+ncor0)] <- init$best[(NEF2+NVC+NW+1):(NEF2+NVC+NW+ncor0)]}
     }
 
     if(ng0==1){
        b <- b1
     }
} else {
    if(length(B)!=NPM){    
       stop("The length of the vector B is not correct")
    } else{
    b <-B
    }
} 

#------------------------------------------
#------nom au vecteur best
#--------------------------------------------



typR <- as.character(typrisq0)
typR[typR=="2"] <- "log(Weibull"
if(logspecif==0){
typR[typR=="1"] <- "+/-sqrt(piecewise"
typR[typR=="3"] <- "+/-sqrt(splines"
}
if(logspecif==1){
typR[typR=="1"] <- "log(piecewise"
typR[typR=="3"] <- "log(splines"
}


if(ng0>=2){
nom <-rep(c("intercept",X0.names[idprob0==1]),each=ng0-1)
   nom1 <-paste(nom," class",c(1:(ng0-1)),sep="")  
   names(b)[1:NPROB]<-nom1
}

hazard.name <- 
if(ng0==1) {
#    names(b)[(1):(nrisqtot)] <- paste("Survparm",c(1:(nprisq0)),sep="")

     names(b)[(1):(nrisqtot)] <- paste(typR,c(1:(nprisq0)),")",sep="")
    names(b)[(nrisqtot+1):(nprisq0+nvarxevt-nvdepsurv)] <- X0.names[idxevt!=0]
	if (nvdepsurv==1)names(b)[(nprisq0+nvarxevt)] <- paste("I(t>",TimeDepVar.name,")",sep="")
    names(b)[(nrisqtot+NVARXEVTinit+1):(NEF)] <- X0.names[idg0!=0]
}


if(ng0>1){
# survie
   if (risqcom0==1){
       names(b)[(NPROB+1):(NPROB+nrisqtot)] <- paste(typR,c(1:(nrisqtot)),")",sep="")
   }


   if (risqcom0==2) {
      names(b)[(NPROB+1):(NPROB+nprisq0)] <- paste(typR,c(1:(nprisq0)),")",sep="")
      names(b)[(NPROB+nprisq0+1):(NPROB+nrisqtot)] <- paste("SurvPH class",c(1:(ng0-1)),sep="")
   }
   if (risqcom0==0){
      for (g in 1:ng0){
	names(b)[((NPROB+nprisq0*(g-1))+1):(NPROB+nprisq0*g)] <- paste(typR,c(1:(nprisq0)),") class",g,sep="")
      }
   }
#noms evol
   nom1<- NULL
   nom2<- NULL
   for (i in 1:nvar.exp) {
      if(idg0[i]==2){ 
          nom <- paste(X0.names[i]," class",c(1:ng0),sep="")	  
          nom1 <- cbind(nom1,t(nom))
      }
      if(idg0[i]==1){
         nom1 <- cbind(nom1,X0.names[i])
      }
# noms survie  
     
#idxevt0
      if(idxevt[i]==2){ 
          nom <- paste(X0.names[i]," class",c(1:ng0),sep="")
          nom2 <- cbind(nom2,t(nom))
      }
      if(idxevt[i]==1) nom2 <- cbind(nom2,X0.names[i])
   }  
    names(b)[(nrisqtot+NPROB+1):(nrisqtot+nvarxevt-nvdepsurv+NPROB)]<- nom2
    if (nvdepsurv==1)names(b)[(nrisqtot+NPROB+nvarxevt)] <- paste("I(t>",TimeDepVar.name,")",sep="")
    names(b)[(nrisqtot+nvarxevt+NPROB+1):NEF]<- nom1

}


if(NVC!=0)names(b)[(NEF+1):(NEF+NVC)] <- paste("Varcov",c(1:(NVC)),sep="")
if(NW!=0)names(b)[(NEF+NVC+1):(NEF+NVC+NW)] <- paste("Varprop class",c(1:(ng0-1)),sep="")
if (ncor0>0) {names(b)[(NEF+NVC+NW+1):(NEF+NVC+NW+ncor0)] <- paste("cor",1:ncor0,sep="") }
names(b)[(NEF+NVC+NW+ncor0+1):(NPM)] <- "Stderr"


N <- NULL
N[1] <- NPROB
N[2] <- nrisqtot
N[3] <- nvarxevt
N[4] <- NEF
N[5] <- NVC
N[6] <- NW
N[7] <- ns1
############################# AJOUT 10/03/2011
N[8] <- nz0
N[9]<- nvdepsurv
N[10] <- ncor0
#############################
idiag <- as.integer(idiag0)
idea <- as.integer(idea0)
nv <- as.integer(nv0)


################ Sortie ###########################


################## Verification des parametres d'entree du modele ##################


 ptm<-proc.time()
cat("Be patient, Jointlcmm is running ... \n")
out <- .Fortran("Jointhet",as.double(Y0),as.double(X0),as.integer(prior0),as.integer(idprob0),as.integer(idea0),as.integer(idg0),as.integer(idcor0),as.integer(idxevt),as.integer(ns0),as.integer(ng0),as.integer(nv0),as.integer(nobs0),as.integer(nea0),as.integer(ncor0),as.integer(nmes0),as.integer(idiag0),as.integer(nwg0),as.integer(NPM),best=as.double(b),V=as.double(V),loglik=as.double(loglik),niter=as.integer(ni),conv=as.integer(istop),gconv=as.double(gconv),ppi2=as.double(ppi0),ppitest2=as.double(ppitest0),resid_m=as.double(resid_m),resid_ss=as.double(resid_ss),pred_m_g=as.double(pred_m_g),pred_ss_g=as.double(pred_ss_g),predRE=as.double(predRE),as.double(convB),as.double(convL),as.double(convG),as.integer(maxiter),as.integer(typrisq0),as.integer(idtrunc0),as.integer(risqcom0),as.integer(nz0),as.double(zi0),as.double(tsurv0),as.double(tsurv),as.double(tsurvint),as.integer(devt),as.integer(indsurvint),statsc=as.double(statsc),risq_est=as.double(risq_est),surv_est=as.double(surv_est),nsim=as.integer(nsim),time_est=as.double(time_est),as.integer(logspecif),PACKAGE="lcmm")


if (!(out$conv %in% c(1,2,4))){
	cat("Problem in the loglikelihood computation. The program stopped abnormally. Please check that the model specification is correct, or try other more plausible initial values.\n")

	out$ppi2 <- rep(NA,ng0*ns0)
	out$ppitest2 <- rep(NA,ng0*ns0)
	out$resid_m <- rep(NA,nobs0)
	out$resid_ss <- rep(NA,nobs0)
	out$pred_m_g <- rep(NA,nobs0*ng0)
	out$pred_ss_g <- rep(NA,nobs0*ng0)
 	out$pred_RE <- rep(NA,ns0*nea0)	
	out$statsc <- NA
	out$risq_est <- rep(NA,nsim*ng0)
	out$surv_est <- rep(NA,nsim*ng0)
	out$V <- rep(NA,(NPM*(NPM+1)/2))
	out$gconv <- rep(NA,3)	

}
if (out$statsc==9999) out$statsc <- NA 

### Creation du vecteur cholesky



Cholesky <- rep(NA,(nea0*(nea0+1)/2))
if(idiag0==0&NVC>0){
   Cholesky[1:NVC] <- out$best[(NEF+1):(NEF+NVC)]
### Construction de la matrice U 
   U <- matrix(0,nrow=nea0,ncol=nea0)
   U[upper.tri(U,diag=TRUE)] <- Cholesky[1:NVC]
   z <- t(U) %*% U
   out$best[(NEF+1):(NEF+NVC)] <- z[upper.tri(z,diag=TRUE)]
}
if(idiag0==1&NVC>0){
   id <- 1:nea0
   indice <- rep(id+id*(id-1)/2)
   Cholesky[indice] <- out$best[(NEF+1):(NEF+nea0)]
   out$best[(NEF+1):(NEF+NVC)] <- out$best[(NEF+1):(NEF+NVC)]**2 
} 


####################################################

if (nea0>0) {
predRE <- matrix(out$predRE,ncol=nea0,byrow=T)
predRE <- data.frame(INDuniq,predRE)
colnames(predRE) <- c(subject.name,X0.names[idea0!=0])
}


if(ng0>1) {
ppi<- matrix(out$ppi2,ncol=ng0,byrow=TRUE)
ppitest<- matrix(out$ppitest2,ncol=ng0,byrow=TRUE)
}
else {
ppi <- matrix(rep(1,ns0),ncol=ng0)
ppitest <- matrix(rep(1,ns0),ncol=ng0)
}


if(!(out$conv %in% c(1,2))){
	classif <-rep(NA,ns0)
}else{
	classif<-apply(ppi,1,which.max)
}


ppi<-data.frame(INDuniq,classif,ppi)
temp<-paste("probYT",1:ng0,sep="")
colnames(ppi) <- c(subject.name,"class",temp)
rownames(ppi) <- 1:ns0
if(!(out$conv %in% c(1,2))){
	classif <-rep(NA,ns0)
}else{
	classif<-apply(ppitest,1,which.max)
}

ppitest<-data.frame(INDuniq,classif,ppitest)
temp<-paste("probY",1:ng0,sep="")
colnames(ppitest) <- c(subject.name,"class",temp)
rownames(ppitest) <- 1:ns0


pred_m_g <- matrix(out$pred_m_g,nrow=nobs0)
pred_ss_g <- matrix(out$pred_ss_g,nrow=nobs0)

if((out$conv %in% c(1,2))){
	pred_m <- Y0-out$resid_m
	pred_ss <- Y0-out$resid_ss
}else{
	pred_m <- rep(NA,nobs0)
	pred_ss <- rep(NA,nobs0)
}


pred <- data.frame(IND,pred_m,out$resid_m,pred_ss,out$resid_ss,Y0,pred_m_g,pred_ss_g)


temp<-paste("pred_m",1:ng0,sep="")
temp1<-paste("pred_ss",1:ng0,sep="")
colnames(pred)<-c(subject.name,"pred_m","resid_m","pred_ss","resid_ss","obs",temp,temp1) 

names(out$best)<-names(b)




surv_est <- matrix(out$surv_est,nrow=nsim)
risq_est <- matrix(out$risq_est,nrow=nsim)
predSurv <- data.frame(out$time_est,risq_est,surv_est)

temp<-paste("RiskFct",1:ng0,sep="")
temp1<-paste("CumRiskFct",1:ng0,sep="")
colnames(predSurv)<-c("time",temp,temp1) 


# rajouter des arguments en sortie
############################# AJOUT survnodes le 10/03/2011

specif <- list(N,ns0,ng0,idea0,idprob0,idg0,idxevt,idiag0,logspecif,idcor0)


survtemp <- update(survival,"~1")
T0.names <- all.names(survtemp)[-c(1,2)]

Names <- list(Y.name,X0.names,T0.names,prior.name,subject.name,X0.names[idea0!=0],TimeDepVar.name)
### ad
Names2 <- list(X0.names2,inddepvar.random.nom2)
### end ad
## ad
res <-list(loglik=out$loglik,best=out$best,V=out$V,gconv=out$gconv,conv=out$conv,call=cl,niter=out$niter,dataset=args$data,pred=pred,pprob=ppi,pprobY=ppitest,predRE=predRE,cholesky=Cholesky,CIstat=out$statsc,predSurv=predSurv,hazard=list(typrisq0,hazardtype,zi0),specif=specif,Names=Names,Names2=Names2,na.action=na.action)

###########################
class(res) <-c("Jointlcmm")
cost<-proc.time()-ptm
cat("The program took", round(cost[3],2), "seconds \n")
res
}
