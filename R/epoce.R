
############################### Fonction Jointlcmm ###################################

epoce <- function(model,pred.times,var.time,newdata=NULL,subset=NULL,na.action=1){

cl <- match.call()
if(missing(var.time)) stop("The argument var.time should be specified")
if (!inherits(var.time, "character")) stop("the class of var.time should be character")
if(missing(model)) stop("The argument model must be specified")
if(class(model)!="Jointlcmm") stop("The argument model must be a class 'Jointlcmm'")


if(!missing(newdata)){
	if(class(newdata)!="data.frame")stop("The argument newdata must be a 'data.frame'")
	if(!all(model$Names[[1]] %in% names(newdata)))stop("The new dataset should contain the dependent longitudinal variable")
### ad
	if(!all(model$Names2[[1]] %in% names(newdata)))stop("The new dataset should contain the all the covariates")
### end ad

	if(!all(model$Names[[3]] %in% names(newdata)))stop("The new dataset should contain the survival variables")
	if(!(var.time %in% names(newdata)))stop("The new dataset should contain the var.time argument")	
	if(!all(model$Names[[5]]  %in% names(newdata)))stop("The new dataset should contain the subject identification variable")	
	if((!(is.null(model$Names[[4]])))&(!all(model$Names[[4]] %in% names(newdata))))stop("The new dataset should contain the Prior variable")		
	new.data <- TRUE
#### Priorname dedans si Priorname diff 0 + Priorname de type integer
} 


if(!(na.action%in%c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument") 

if(na.action==1){
	na.action=na.omit
}else{
	na.action=na.fail
}



######### RECUP SPECIFICATION ##########################
nT <- length(pred.times)
rl_cond <- rep(0,nT)
epoir <- rep(0,nT)
vopt <- as.double(model$V)
best <- as.double(model$best)
NPM <- length(best)
ns_vect <- rep(0,nT)
nevt_vect <- rep(0,nT)

### specification recup de args
idprob0 <- model$specif[[5]]
idea0 <- model$specif[[4]]
idg0 <- model$specif[[6]]
idcor0 <- model$specif[[10]]
idxevt <- model$specif[[7]]
idxevt <- model$specif[[7]]
idiag0 <- model$specif[[8]]
nv0 <- length(idprob0)
nwg0 <- model$specif[[1]][6]
ng0 <- model$specif[[3]]
ncor0 <- model$specif[[1]][10]
nz0 <- model$specif[[1]][8]
zi0 <- model$hazard[[3]]
typrisq0 <-  model$hazard[[1]]
risqcom0 <- switch(model$hazard[[2]],"Specific"=0,"PH"=2,"Common"=1)
logspecif <- model$specif[[9]]

#cat(paste("  risqcom ",risqcom0)," \n")
#cat(paste("  idprob0 ",idprob0)," \n")
#cat(paste("  idea0 ",idea0)," \n")
#cat(paste("  idg0 ",idg0)," \n")
#cat(paste("  idpxevt ",idxevt)," \n")

##############   RECUP DATA   ##########################


### new data or data from estimation

if(missing(newdata)){
	data <- eval(model$data)
	if(!(var.time %in% names(data)))stop("The Jointlcmm data must be contain the var.time variable")
	new.data <- FALSE
###  var.time est dedans
}else{
	data <- newdata
}

int.fixed <- 0
int.mixture <- 0
int.random <- 0
int.classmb <- 0
#7/05/2012
### Traitement des donnees manquantes
# fixed
mcall <- model$call[c(1,match(c("data"),names(model$call),0))]
mcall$na.action <- na.action
mcall$subset <- subset
mcall$data <- data	

m <- mcall
m$formula <- model$call$fixed
m[[1]] <- as.name("model.frame")	
m <- eval(m, sys.parent()) 
na.fixed <- attr(m,"na.action")

# mixture
if(!is.null(model$call$mixture)){
	m <- mcall
	m$formula <- model$call$mixture
	m[[1]] <- as.name("model.frame")	
	m <- eval(m, sys.parent()) 
	na.mixture <- attr(m,"na.action")	
}else{
	na.mixture <- NULL
}
# random
if(!is.null(model$call$random)){
	m <- mcall
	m$formula <- model$call$random
	m[[1]] <- as.name("model.frame")	
	m <- eval(m, sys.parent()) 
 	na.random <- attr(m,"na.action")
}else{
	na.random <- NULL
}
# classmb
if(!is.null(model$call$classmb)){ 
	m <- mcall	
	m$formula <- model$call$classmb	
	m[[1]] <- as.name("model.frame")	
	m <- eval(m, sys.parent()) 
 	na.classmb <- attr(m,"na.action")
}else{
	na.classmb <- NULL
}
#7/05/2012
########## For survival
if(!is.null(model$call$survival)){ 
	res.evt <- model$call$survival
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
	m <- mcall
	m$formula <- update(tmp.res,names.survival)
	m[[1]] <- as.name("model.frame")
	m <- eval(m, sys.parent())
	na.survival <- attr(m,"na.action")

## Table sans donnees manquante: newdata
	na.action <- unique(c(na.fixed,na.mixture,na.random,na.classmb,na.survival))
	if(!is.null(na.action)){
		newdata1 <- data[-na.action,]
	}else{
		newdata1 <- data
	}

#7/05/2012
	X_survival <- model.matrix(formula(update(tmp.res,names.survival)),data=newdata1)

	if(colnames(X_survival)[1]=="(Intercept)"){
		colnames(X_survival)[1] <- "intercept"
		int.survival <- 1
	}
	if(length(inddepvar.Mixt) > 0){
		inddepvar.Mixt <- colnames(X_survival)[grep(inddepvar.Mixt,colnames(X_survival))]
		inddepvar.noMixt <- colnames(X_survival)[!(colnames(X_survival) %in% inddepvar.Mixt)][-1]
	}else{
		inddepvar.noMixt <- colnames(X_survival)[-1]
	}
	if(int.survival>0) inddepvar.survival <- colnames(X_survival)[-1]
	id.X_survival <- 1
}else{
	id.X_survival <- 0
#7/05/2012
	na.action <- unique(c(na.fixed,na.mixture,na.random,na.classmb))
	if(!is.null(na.action)){
		newdata1 <- data[-na.action,]
	}else{
		newdata1 <- data
	}
#7/05/2012
}

### DEP VAR
Y0 <- newdata1[,model$Names[[1]]]
nobs0 <- length(Y0)
### SURVIVAL TIMES

if (length(model$Names[[3]])==2) {
Tevent <- newdata1[,model$Names[[3]][1]]
Devent <- newdata1[,model$Names[[3]][2]]
Tentry <- rep(0,nobs0)
}
if (length(model$Names[[3]])==3) {
Tentry <- newdata1[,model$Names[[3]][1]]
Tevent <- newdata1[,model$Names[[3]][2]]
Devent <- newdata1[,model$Names[[3]][3]]
}

if ((max(Tevent)>max(model$hazard[[3]]))&(model$hazard[[1]]!=2)) stop("The maximal time of event in the new dataset should not be greater than the maximal time of event in the dataset used in Jointlcmm.")

#### TimeDepVar
nvdepsurv <- model$specif[[1]][9]
Tint <- Tevent
indsurvint <- rep(0,length= model$specif[[2]])
if(nvdepsurv!=0)  {
Tint <- newdata1[,model$Names[[7]]]
Tint[(is.na(Tint))] <- Tevent[(is.na(Tint))]
Tint[Tint>Tevent] <- Tevent[Tint>Tevent]
Tint[Tint<Tentry] <- Tentry[Tint<Tentry]
if (length(Tint[Tint>Tentry])==0) Tint <- Tevent
}
## Construction de nouvelles var eplicatives sur la nouvelle table
## fixed
	
	X_fixed <- model.matrix(formula(model$call$fixed),data=newdata1)
	if(colnames(X_fixed)[1]=="(Intercept)"){
		colnames(X_fixed)[1] <- "intercept"
		int.fixed <- 1
	}
## mixture
	if(!is.null(model$call$mixture)){
		X_mixture <- model.matrix(formula(model$call$mixture),data=newdata1)	
		if(colnames(X_mixture)[1]=="(Intercept)"){
			colnames(X_mixture)[1] <- "intercept"
			int.mixture <- 1
		}
		id.X_mixture <- 1
	}else{
		id.X_mixture <- 0
	}
## random
	if(!is.null(model$call$random)){
		X_random <- model.matrix(formula(model$call$random),data=newdata1)	
		if(colnames(X_random)[1]=="(Intercept)"){
			colnames(X_random)[1] <- "intercept"
			int.random <- 1
		}
		id.X_random <- 1
	}else{
		id.X_random <- 0
	}
## classmb
	if(!is.null(model$call$classmb)){ 
		X_classmb <- model.matrix(formula(model$call$classmb),data=newdata1)
		colnames(X_classmb)[1] <- "intercept"
		id.X_classmb <- 1
	}else{
		id.X_classmb <- 0
	}

varX0.names <- NULL
varX0.names <- c(varX0.names,colnames(X_fixed))
if(id.X_mixture == 1) varX0.names <- c(varX0.names,colnames(X_mixture))
if(id.X_random == 1) varX0.names <- c(varX0.names,colnames(X_random))
if(id.X_classmb == 1) varX0.names <- c(varX0.names,colnames(X_classmb))
varX0.names <- unique(varX0.names)


if(id.X_survival == 1) varX0.names <- unique(c(varX0.names,colnames(X_survival)))
### construction de X0
## var expli
X0 <- X_fixed
if(id.X_mixture == 1){
	for(i in 1:length(colnames(X_mixture))){
		if((colnames(X_mixture)[i] %in% colnames(X0))==F){
			X0 <- cbind(X0,X_mixture[,i])
			
		}
	}
}
if(id.X_random == 1){
	for(i in 1:length(colnames(X_random))){
		if((colnames(X_random)[i] %in% colnames(X0))==F){
			X0 <- cbind(X0,X_random[,i])
		}	 
	}
}
if(id.X_classmb == 1){
	for(i in 1:length(colnames(X_classmb))){
		if((colnames(X_classmb)[i] %in% colnames(X0))==F){
			X0 <- cbind(X0,X_classmb[,i],deparse.level=0)	 
		}	
	}
}
if(id.X_survival == 1){
	for(i in 1:length(colnames(X_survival))){
		if((colnames(X_survival)[i] %in% colnames(X0))==F){
			X0 <- cbind(X0,X_survival[,i])
			
		}
	}
}
colnames(X0) <- varX0.names
X0 <- as.data.frame(X0)


### end ad 2/04/2012

### VAR EXP : attention intercept or no 
# if (!("intercept" %in% model$Names[[2]])){
# 	X0 <- as.data.frame(data[,model$Names[[2]]]) 
# }
# if ("intercept" %in% model$Names[[2]]){
# 	Xnames2 <- setdiff(model$Names[[2]],"intercept")
# 	X0 <- as.data.frame(data[,Xnames2]) 
# 	X0<- cbind(intercept=rep(1,nobs0),X0)
# 	names(X0) <- model$Names[[2]]
# }
### end ad 2/04/2012

### IND

if(missing(newdata)){
	IND <- newdata1[,colnames(model$pprob)[1]]
}else{
	IND <- newdata1[,model$Names[[5]]]
}
IDnum <- as.numeric(IND)

### Time
Time<-newdata1[,var.time]


#### PRIOR : A REVOIR EN ENTIER ?????????????????

##### INCLUSION PRIOR 
if(is.null(model$Priorname)){ 
PRIOR <- seq(0,length=length(IND))} 
else
{
PRIOR <- newdata1[,model$Priorname]
PRIOR[(is.na(PRIOR))] <- 0
}


###### DATA SORTING on IND variable
matYX <- cbind(IDnum,IND,PRIOR,Tentry,Tevent,Devent,Tint,Y0,Time,X0)
matYXord <- matYX[sort.list(matYX[,1]),]
Y0 <- matYXord[,8]
Time <- matYXord[,9]
X0 <- matYXord[,-c(1,2,3,4,5,6,7,8,9)]
IDnum <- matYXord[,1]
IND <- matYXord[,2]
PRIOR <- matYXord[,3]
PRIOR <- as.integer(as.vector(PRIOR))
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

INDuniq <- IND[cumsum(nmes0)]

#cat(paste("vopt",length(vopt))," \n")
#cat(paste(vopt),"\n")

#cat(paste("Best",length(best))," \n")
#cat(paste(best),"\n")

#cat(paste("pred.times",length(pred.times))," \n")
#cat(paste(pred.times),"\n")

#cat(paste("nmes0",length(nmes0))," \n")
#cat(paste(nmes0)," \n")
#cat("head matYX \n")
#cat(paste(matYXord[1,])," \n")
#cat("head(X0) \n")
#cat("avant prior \n")

### reduction de prior
initial <- as.integer(rep(0,ns0))
prior0 <- initial
# si prior pas missing alors mettre dedans la classe a priori. Attention tester q les valeurs sont dans 0, G
#prior0 <- PRIOR[cumsum(nmes0)]


initial <- as.integer(rep(0,ns0))

### reduction de Devent, Tevent,Tentry a la taille ns0
devt <- initial
tsurv0 <- initial
tsurv <- initial
tsurvint <- initial
devt <- Devent[cumsum(nmes0)]
tsurv <- Tevent[cumsum(nmes0)]
tsurv0 <- Tentry[cumsum(nmes0)]
tsurvint <- Tint[cumsum(nmes0)]
indsurvint[tsurvint<tsurv] <- 1

idtrunc0 <- 0
if (all.equal(tsurv0,0)==F) {idtrunc0 <-1}

#cat(paste("devt",c(length(devt),sum(devt)/length(devt)))," \n")
#cat(paste("tsurv",c(min(tsurv),max(tsurv)))," \n")
#cat(paste("tsurv",c(min(tsurv0),max(tsurv0)))," \n")

contribt <- rep(0,length=ns0*nT)

################ FORTRAN FUNCTION CALL #####################
#
ptm<-proc.time()
cat("Be patient, epoce function is running ... \n")
 
#cat("c(nT,ns0,nobs0,ng0,nv0,idiag0,nwg0,NPM,typrisq0,idtrunc0,risqcom0,nz0,nvdepsurv)"," \n")
#cat(paste(c(nT,ns0,nobs0,ng0,nv0,idiag0,nwg0,NPM,typrisq0,idtrunc0,risqcom0,nz0,nvdepsurv))," \n")


#cat(paste("zi0",zi0)," \n")
#cat(paste("ns_vect",ns_vect)," \n")
#cat(paste("nevt_vect",nevt_vect)," \n")

#cat(paste("specif",logspecif)," \n")


out <- .Fortran("cvpl",as.double(Y0),as.double(X0),as.integer(prior0),as.integer(idprob0),as.integer(idea0),as.integer(idg0),as.integer(idcor0),as.integer(idxevt),as.integer(ns0),as.integer(ng0),as.integer(ncor0),as.integer(nv0),as.integer(nobs0),as.integer(nmes0),as.integer(idiag0),as.integer(nwg0),as.integer(NPM),as.double(Time),as.integer(typrisq0),as.integer(idtrunc0),as.integer(risqcom0),as.integer(nz0),as.double(zi0),as.double(tsurv0),as.double(tsurv),as.double(tsurvint),as.integer(devt),as.integer(indsurvint),as.double(vopt),as.integer(nT),as.double(pred.times),best=as.double(best),epoir=as.double(epoir),rl_cond=as.double(rl_cond),ns_vect=as.integer(ns_vect),nevt_vect=as.integer(nevt_vect),contribt=as.double(contribt),as.integer(logspecif),PACKAGE="lcmm")

# construction de la matrice contribt

contrib <- matrix(out$contribt,nrow=ns0,ncol=nT)
namesContrib <- as.vector(apply(matrix(pred.times,nrow=1),MARGIN=2,FUN=function(x){paste("IndivContrib_time_",x,sep="")}))
colnames(contrib) <- namesContrib
contrib <- cbind(INDuniq,contrib)
contrib[contrib==0] <- NA
# remplacer les 0 (vrais 0 par NA) (le nombre de non nul pour un temps = ns_vect(de ce temps)
# colnames : permier colonne = le nom de la variable IND = colnames(model$pprob)[1]
#          : colonnes suivantes = IndivContrib_time_x, x etant la valeur du temps de prediction (ce qu'il y a dans pred.times)

if (!is.null(newdata)){
out$epoir <- rep(NA,length(pred.times))
}
out$epoir[out$epoir==1.e9] <- NA
out$rl[out$rl==-1.e9] <- NA


cvpl <- cbind(pred.times,out$ns_vect,out$nevt_vect,-out$rl,out$epoir)
colnames(cvpl) <- c("pred. times"," N at risk","N events","MPOL","CVPOL")
rownames(cvpl) <- rep(" ",length(cvpl[,1]))

# sortie des resultats
res <- list(call.Jointlcmm=model$call,call.epoce=cl,EPOCE=cvpl,IndivContrib=contrib,new.data=new.data)

class(res) <-c("epoce")
cost<-proc.time()-ptm
cat("The program took", round(cost[3],2), "seconds \n")
res
}
