###################### derniere mise a jour : 2012/03/16 ############"
.Ordlcmm <-
function(fixed,mixture,random,subject,classmb,ng,idiag,nwg,data,B,convB,convL,convG,prior,maxiter,zitr,ide,call,Ydiscrete,subset=subset,na.action)

{

ptm<-proc.time()
cat("Be patient, lcmm is running ... \n")

cl <- match.call()

args <- as.list(match.call(.Ordlcmm))[-1]

#nom.subject <- as.character(args$subject)
nom.subject <- as.character(subject)

#### INCLUSION PRIOR
nom.prior <- as.character(args$prior)
####
if(!missing(mixture) & ng==1) stop("No mixture can be specified with ng=1")
if(missing(mixture) & ng>1) stop("The argument mixture has to be specified for ng > 1")
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

if(missing(data)){ stop("The argument data should be specified and defined as a data.frame")} 
if(missing(subject)){ stop("The argument subject must be specified in any model even without random-effects")} 
### ad 2/04/2012
X0.names2 <- c("intercept")
### ad
int.fixed <- 0
int.mixture <- 0
int.random <- 0
int.classmb <- 0

## Table sans donnees manquante: newdata
newdata <- data  
#prendre le subset :
if(!(is.null(subset))) 
newdata <- data[subset,]

#enlever les NA
	if(!is.null(na.action)){
		newdata <- newdata[-na.action,]
	}

### names of covariate in intial fit
X0.names2 <- unique(c(X0.names2,colnames(get_all_vars(formula(terms(fixed)),data=newdata))[-1]))
if(mixture[[2]] != "-1")X0.names2 <- unique(c(X0.names2,colnames(get_all_vars(formula(terms(mixture)),data=newdata))))
if(random[[2]] != "-1")X0.names2 <- unique(c(X0.names2,colnames(get_all_vars(formula(terms(random)),data=newdata))))
if(classmb[[2]] != "-1")X0.names2 <- unique(c(X0.names2,colnames(get_all_vars(formula(terms(classmb)),data=newdata))))


## Construction de nouvelles var eplicatives sur la nouvelle table
## fixed
	X_fixed <- model.matrix(fixed,data=newdata)
	if(colnames(X_fixed)[1]=="(Intercept)"){
		colnames(X_fixed)[1] <- "intercept"
		int.fixed <- 1
	}else{
		stop ("Only models with an intercept can be estimated using lcmm. This is required for identifiability purposes")
	}
	nom.fixed <- colnames(X_fixed)
	if(int.fixed>0)inddepvar.fixed <- inddepvar.fixed.nom <- nom.fixed[-1]

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
		inddepvar.random <- inddepvar.random.nom <- NULL
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

# intercept is always in inddepvar.classmb
var.exp <- NULL
var.exp <- c(var.exp,colnames(X_fixed))
if(id.X_mixture == 1) var.exp <- c(var.exp,colnames(X_mixture))
if(id.X_random == 1)var.exp <- c(var.exp,colnames(X_random))
if(id.X_classmb == 1)var.exp <- c(var.exp,colnames(X_classmb))
var.exp <- unique(var.exp)
## ad

#if(!(all(nom.mixture %in% nom.fixed))) stop("The covariates in mixture should be also included in the argument fixed")
# controler si les variables de mixture sont toutes dans fixed : 
 z.fixed <- strsplit(nom.fixed,split=":",fixed=TRUE)
 z.fixed <- lapply(z.fixed,sort)
  
 if(id.X_mixture==1)
 {
  z.mixture <- strsplit(nom.mixture,split=":",fixed=TRUE)
  z.mixture <- lapply(z.mixture,sort)
 }
 else z.mixture <- list()

if(!all(z.mixture %in% z.fixed))  stop("The covariates in mixture should also be included in the argument fixed")

#ad 
## var dependante
Y.name <- as.character(attributes(terms(fixed))$variables[2])
Y0 <- newdata[,Y.name]
#ad 

## var expli
X0 <- X_fixed
oldnames <- colnames(X0)

z.X0 <- strsplit(colnames(X0),split=":",fixed=TRUE)
z.X0 <- lapply(z.X0,sort)

if(id.X_mixture == 1)
{
 z.mixture <- strsplit(colnames(X_mixture),split=":",fixed=TRUE)
 z.mixture <- lapply(z.mixture,sort)
	for(i in 1:length(colnames(X_mixture)))
 {
		#if((colnames(X_mixture)[i] %in% colnames(X0))==F)
  if(!isTRUE(z.mixture[i] %in% z.X0))
  {
			X0 <- cbind(X0,X_mixture[,i])
			colnames(X0) <- c(oldnames, colnames(X_mixture)[i])
			oldnames <- colnames(X0)
   			
   z.X0 <- strsplit(colnames(X0),split=":",fixed=TRUE)
   z.X0 <- lapply(z.X0,sort)					
		}                                                                
	}
}
else
{
 z.mixture <- list()
}

if(id.X_random == 1)
{
 z.random <- strsplit(colnames(X_random),split=":",fixed=TRUE)
 z.random <- lapply(z.random,sort)
	for(i in 1:length(colnames(X_random)))
 {
		#if((colnames(X_random)[i] %in% colnames(X0))==F)
		if(!isTRUE(z.random[i] %in% z.X0))
  {
			X0 <- cbind(X0,X_random[,i])
			colnames(X0) <- c(oldnames, colnames(X_random)[i])
			oldnames <- colnames(X0)
   			
   z.X0 <- strsplit(colnames(X0),split=":",fixed=TRUE)
   z.X0 <- lapply(z.X0,sort)			
		}	 
	}
}
else
{
 z.random <- list()
}

if(id.X_classmb == 1)
{
 z.classmb <- strsplit(colnames(X_classmb),split=":",fixed=TRUE)
 z.classmb <- lapply(z.classmb,sort)
	for(i in 1:length(colnames(X_classmb)))
 {
		#if((colnames(X_classmb)[i] %in% colnames(X0))==F)
		if(!isTRUE(z.classmb[i] %in% z.X0))
  {
			X0 <- cbind(X0,X_classmb[,i])
			colnames(X0) <- c(oldnames, colnames(X_classmb)[i])
			oldnames <- colnames(X0)
      	
   z.X0 <- strsplit(colnames(X0),split=":",fixed=TRUE)
   z.X0 <- lapply(z.X0,sort)    
		}	
	}
}
else
{
 z.classmb <- list()
}

#colnames(X0) <- var.exp


if((any(is.na(X0))==TRUE)|(any(is.na(Y0))==TRUE))stop("The data should not contain any missing value")
 
n <- dim(data)[1]
#if ((int.fixed+int.random)>0) X0<- cbind(intercept=rep(1,n),X0)
if (!((int.fixed+int.random)>0)) X0 <- as.data.frame(X0[,-which(colnames(X0)=="intercept")])
nom.X0 <- colnames(X0)
nvar.exp <- length(nom.X0)

IND <- newdata[,nom.subject]
IDnum <- as.numeric(IND)


#### INCLUSION PRIOR 
if(missing(prior)){ PRIOR <- seq(0,length=length(IND))} 
if(!missing(prior)){ 
PRIOR <- newdata[,nom.prior]
PRIOR[(is.na(PRIOR))] <- 0
}
####

ng0 <- ng
idiag0 <- as.integer(idiag)
nwg0 <- as.integer(nwg)

idea0 <- rep(0,nvar.exp)
idprob0 <- rep(0,nvar.exp)
idg0 <- rep(0,nvar.exp)

z.X0 <- strsplit(nom.X0,split=":",fixed=TRUE)
z.X0 <- lapply(z.X0,sort)

for (i in 1:nvar.exp)    
{
# idea0[i] <- nom.X0[i]%in%inddepvar.random.nom
# idprob0[i] <- nom.X0[i]%in%inddepvar.classmb.nom      
# if(nom.X0[i]%in%nom.fixed & !(nom.X0[i]%in%nom.mixture)) idg0[i] <- 1 
# if(nom.X0[i]%in%nom.fixed & nom.X0[i]%in%nom.mixture) idg0[i] <- 2 
 
 idea0[i] <- z.X0[i] %in% z.random
 idprob0[i] <- z.X0[i] %in% z.classmb   
 if((z.X0[i] %in% z.fixed) & !(z.X0[i] %in% z.mixture)) idg0[i] <- 1 
 if((z.X0[i] %in% z.fixed) & (z.X0[i] %in% z.mixture)) idg0[i] <- 2     
}

if((int.fixed+int.random)>0) idprob0[1] <- 0

# on ordonne les donn es suivants la variable IND
matYX <- cbind(IDnum,IND,PRIOR,Y0,X0)
matYXord <- matYX[sort.list(matYX[,1]),]
Y0 <- matYXord[,4]  
X0 <- matYXord[,-c(1,2,3,4)]
IDnum <- matYXord[,1]
IND <- matYXord[,2]


#### INCLUSION PRIOR 
PRIOR <- matYXord[,3]
PRIOR <-as.integer(as.vector(PRIOR))
####

X0<-as.numeric(as.matrix(X0))
Y0<-as.numeric(as.matrix(Y0))
nmes0<-as.vector(table(IND))
ns0<-length(nmes0)


##### INCLUSION PRIOR 
# definition de prior a 0 pour l'analyse G=1
prior2 <- as.integer(rep(0,ns0))
prior0 <- prior2 
# si prior pas missing alors mettre dedans la classe a priori. Attention tester q les valeurs sont dans 0, G
if(!missing(prior)){ 
prior0 <- PRIOR[cumsum(nmes0)]
}
INDuniq <- IND[cumsum(nmes0)]
seqnG <- 0:ng0
if (!(all(prior0  %in% seqnG))) stop ("The argument prior should contain integers between 0 and ng")
#####


loglik <- as.double(0)
rlindiv <- rep(0,ns0)
UACV <- as.double(0)
ni <- 0
istop <- 0
gconv <-rep(0,3)
ppi0 <- rep(0,ns0*ng0)
nv0<-nvar.exp
nobs0<-length(Y0)
resid_m <- rep(0,nobs0)
resid_ss <- rep(0,nobs0)
pred_m_g <- rep(0,nobs0*ng0)
pred_ss_g <- rep(0,nobs0*ng0)
nea0 <- sum(idea0==1)
predRE <- rep(0,nea0*ns0)

ntrtot0 <- sum(ide==1) 
minY <- zitr[1]
maxY <- zitr[2]





#-------------------------------------------------------------------------------
#definition du vecteur de parametre + initialisation
#-------------------------------------------------------------------------------
#####cas 1 : ng=1
b<-NULL
b1 <- NULL
NPROB <- 0
if(ng0==1| missing(B)){
NEF<-sum(idg0!=0)-1
b1[1:NEF]<-0

NVC <- 0
if(idiag0==1&nea0>0){
NVC<-sum(idea0==1)
b1[(NEF+1):(NEF+NVC)]<-1}

if(idiag0==0&nea0>0){
kk<-sum(idea0==1) 
NVC<-(kk*(kk+1))/2
indice<-cumsum(1:kk)
bidiag<-rep(0,NVC)
bidiag[indice]<-1
b1[(NEF+1):(NEF+NVC)]<-bidiag
}

# valeurs initiales pour les transfos #
if (ntrtot0==1) {
b1[(NEF+NVC+1):(NEF+NVC+ntrtot0)] <- 0
}else{
b1[(NEF+NVC+1)]<- 2*qnorm(0.98)*(-median(Y0)+minY+1)/(ntrtot0-1)
b1[(NEF+NVC+2):(NEF+NVC+ntrtot0)] <- sqrt(2*qnorm(0.98)/(ntrtot0-1))
}

NPM<-length(b1)
NW<-0
V <- rep(0,NPM*(NPM+1)/2) 
}

#####cas 2 : ng>=2
if(ng0>1){
NPROB<-(sum(idprob0==1)+1)*(ng0-1)
b[1:NPROB]<-0
NEF<-sum(idg0==1)+(sum(idg0==2))*ng0-1
NVC <- 0
if(idiag0==1&nea0>0)NVC<-sum(idea0==1)
if(idiag0==0&nea0>0){
kk<-sum(idea0==1) 
NVC<-(kk*(kk+1))/2}
NW<-nwg0*(ng0-1)
if(NW>0) b[(NPROB+NEF+NVC+1):(NPROB+NEF+NVC+NW)]<-1
NPM<-NPROB+NEF+NVC+NW+ntrtot0
V <- rep(0,NPM*(NPM+1)/2)
} 



if(missing(B)){

if(ng0>1){
idea2 <- idea0
idprob2 <- rep(0,nv0)  
idg2 <- rep(0,nv0) 
idg2[idg0!=0] <- 1
NEF2<-sum(idg2==1)-1
NPM2<-NEF2+NVC+ntrtot0
nwg2<-0
ng2<-1
ppi2<- rep(0,ns0)
pred_m_g2 <- rep(0,nobs0)
pred_ss_g2 <- rep(0,nobs0)

V2 <- rep(0,NPM2*(NPM2+1)/2)

marker <- rep(0,(maxY-minY+1)*2)
transfY <- rep(0,(maxY-minY+1)*2)
init <- .Fortran("hetmixOrd",as.double(Y0),as.double(X0),as.integer(prior2),as.integer(idprob2),
as.integer(idea2),as.integer(idg2),as.integer(ns0),as.integer(ng2),as.integer(nv0),as.integer(nobs0),
as.integer(nea0),as.integer(nmes0),as.integer(idiag0),as.integer(nwg2),as.integer(NPM2),
best=as.double(b1),V=as.double(V2),as.double(loglik),niter=as.integer(ni),conv=as.integer(istop),
as.double(gconv),as.double(ppi2),as.double(resid_m),as.double(resid_ss),as.double(pred_m_g2),
as.double(pred_ss_g2),predRE=as.double(predRE),as.double(convB),as.double(convL),as.double(convG),
as.integer(maxiter),as.integer(minY),as.integer(maxY),as.integer(ide),as.double(marker),as.double(transfY),
as.double(UACV),as.double(rlindiv),PACKAGE="lcmm")


k <- NPROB
l <- 0
t<- 0
for (i in 1:nvar.exp)    {
if(idg0[i]==1 & i>1){
l <- l+1
t <- t+1
b[k+t] <- init$best[l]
}
if(idg0[i]==2){
if (i==1){
	for (g in 2:ng){
	t <- t+1
	b[k+t] <- - 0.5*(g-1)
	}
}
if (i>1){
	l <- l+1
	for (g in 1:ng){
	t <- t+1
	if(init$conv==1) b[k+t] <- init$best[l]+(g-(ng+1)/2)*sqrt(init$V[l*(l+1)/2])
	else  b[k+t] <- init$best[l]+(g-(ng+1)/2)*init$best[l]
	}
}
}
}

b[(NPROB+NEF+1):(NPROB+NEF+NVC)] <-init$best[(NEF2+1):(NEF2+NVC)]
b[(NPM-ntrtot0+1):NPM] <-init$best[(NPM2-ntrtot0+1):NPM2]
} 
if(ng0==1 ){
b <- b1
}
} else {if(length(B)!=NPM){
stop("The length of the vector B is not correct")
}
 else {b <-B}
} 

se <- rep(0,length(b))
#------------------------------------------
#------nom au vecteur best
#--------------------------------------------

#if(ng0==2)names(b)[1:NPROB]<-nom.X0[idprob0!=0]

if(ng0>=2){
nom <-rep(c("intercept",nom.X0[idprob0!=0]),each=ng0-1)
nom1 <- paste(nom," class",c(1:(ng0-1)),sep="")
names(b)[1:NPROB]<-nom1
}



if(ng0==1) names(b)[1:(NEF)] <-  nom.X0[-1][idg0[-1]!=0]

if(ng0>1){
	nom1<- NULL
	for (i in 1:nvar.exp) {
		if(idg0[i]==2){ 
		   if (i==1){
			 nom <- paste(nom.X0[i]," class",c(2:ng0),sep="")
		       nom1 <- cbind(nom1,t(nom))
		    }
		   if (i>1){
			 nom <- paste(nom.X0[i]," class",c(1:ng0),sep="")
			 nom1 <- cbind(nom1,t(nom))
		    }
		}
	      if(idg0[i]==1 & i>1) nom1 <- cbind(nom1,nom.X0[i])
	}
names(b)[(NPROB+1):(NPROB+NEF)]<- nom1
}

names(b)[(NPM-ntrtot0+1):NPM]<- paste("thresh. parm",c(1:ntrtot0),sep="")



if(NVC!=0)names(b)[(NPROB+NEF+1):(NPROB+NEF+NVC)] <- paste("varcov",c(1:(NVC)))
if(NW!=0)names(b)[(NPROB+NEF+NVC+1):(NPROB+NEF+NVC+NW)] <- paste("varprop class",c(1:(ng0-1)))

N <- NULL
N[1] <- NPROB
N[2] <- NEF
N[3] <- NVC
N[4] <- NW
N[5] <- nobs0
N[6] <- 0 #=ncor

idiag <- as.integer(idiag0)
idea <- as.integer(idea0)
nv <- as.integer(nv0)


################ Sortie ###########################


marker <- rep(0,(maxY-minY+1)*2)
transfY <- rep(0,(maxY-minY+1)*2)

#cat("avant le programme\n")



out <- .Fortran("hetmixOrd",as.double(Y0),as.double(X0),as.integer(prior0),as.integer(idprob0),
as.integer(idea0),as.integer(idg0),as.integer(ns0),as.integer(ng0),as.integer(nv0),as.integer(nobs0),
as.integer(nea0),as.integer(nmes0),as.integer(idiag0),as.integer(nwg0),as.integer(NPM),best=as.double(b),
V=as.double(V),loglik=as.double(loglik),niter=as.integer(ni),conv=as.integer(istop),gconv=as.double(gconv),
ppi2=as.double(ppi0),resid_m=as.double(resid_m),resid_ss=as.double(resid_ss),pred_m_g=as.double(pred_m_g),
pred_ss_g=as.double(pred_ss_g),predRE=as.double(predRE),as.double(convB),as.double(convL),as.double(convG),
as.integer(maxiter),as.integer(minY),as.integer(maxY),as.integer(ide),marker=as.double(marker),
transfY=as.double(transfY),UACV=as.double(UACV),rlindiv=as.double(rlindiv),PACKAGE="lcmm")

### Creation du vecteur cholesky


Cholesky <- rep(0,(nea0*(nea0+1)/2))
if(idiag0==0 & NVC>0){
Cholesky[1:NVC] <- out$best[(NPROB+NEF+1):(NPROB+NEF+NVC)]
### Construction de la matrice U 
U <- matrix(0,nrow=nea0,ncol=nea0)
U[upper.tri(U,diag=TRUE)] <- Cholesky[1:NVC]
z <- t(U) %*% U
out$best[(NPROB+NEF+1):(NPROB+NEF+NVC)] <- z[upper.tri(z,diag=TRUE)]
}
if(idiag0==1 & NVC>0){
id <- 1:nea0
indice <- rep(id+id*(id-1)/2)
Cholesky[indice] <- out$best[(NPROB+NEF+1):(NPROB+NEF+nea0)]
out$best[(NPROB+NEF+1):(NPROB+NEF+NVC)] <- out$best[(NPROB+NEF+1):(NPROB+NEF+NVC)]**2 
} 

####################################################

#predRE <- matrix(out$predRE,ncol=nea0,byrow=T)
#predRE <- cbind(INDuniq,predRE)
#colnames(predRE) <- c(nom.subject,nom.X0[idea0!=0])


if(ng0>1) {
ppi<- matrix(out$ppi2,ncol=ng0,byrow=TRUE)
}
else {
ppi <- matrix(rep(1,ns0),ncol=ng0)
}


classif<-apply(ppi,1,which.max)
ppi<-data.frame(INDuniq,classif,ppi)
temp<-paste("prob",1:ng0,sep="")
colnames(ppi) <- c(nom.subject,"class",temp)
rownames(ppi) <- 1:ns0

#pred_m_g <- matrix(out$pred_m_g,nrow=nobs0)
#pred_ss_g <- matrix(out$pred_ss_g,nrow=nobs0)
#pred_m <- Y0-out$resid_m
#pred_ss <- Y0-out$resid_ss
#pred <- cbind(IND,pred_m,out$resid_m,pred_ss,out$resid_ss,Y0,pred_m_g,pred_ss_g)

#temp<-paste("pred_m",1:ng0,sep="")
#temp1<-paste("pred_ss",1:ng0,sep="")
#colnames(pred)<-c(nom.subject,"pred_m","resid_m","pred_ss","resid_ss","obs",temp,temp1) 

names(out$best)<-names(b)
#btest <- out$best[1:length(inddepvar.fixed.nom)]
#names(btest) <-inddepvar.fixed.nom


estimlink <- cbind(out$marker,out$transfY)
colnames(estimlink) <- c("Y","transfY")

### ad 2/04/2012
if (!("intercept" %in% nom.X0)) X0.names2 <- X0.names2[-1]
### ad

res <-list(ns=ns0,ng=ng0,idea0=idea0,idprob0=idprob0,idg0=idg0,idcor0=rep(0,nvar.exp),loglik=out$loglik,best=out$best,V=out$V,gconv=out$gconv,conv=out$conv,call=call,niter=out$niter,N=N,idiag=idiag0,pprob=ppi,Xnames=nom.X0,Xnames2=X0.names2,cholesky=Cholesky,estimlink=estimlink,linktype=3,linknodes=zitr,ide=ide,Ydiscrete=Ydiscrete,discrete_loglik=out$loglik,UACV=out$UACV,IndivContrib=out$rlindiv,na.action=na.action)
class(res) <-c("lcmm")  

cost<-proc.time()-ptm
cat("The program took", round(cost[3],2), "seconds \n")
res
}


