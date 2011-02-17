
############################### Fonction hlme ###################################

hlme <-
function(fixed,mixture,random,subject,classmb,ng=1,idiag=FALSE,nwg=FALSE,data,B,convB=0.0001,convL=0.0001,convG=0.0001,prior,Maxiter=100){
cat("Be patient. The program is running ... \n")

cl <- match.call()
args <- as.list(match.call(hlme))[-1]

nom.subject <- as.character(args$subject)
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

###########################################################"
res.fixed <- terms(fixed)
attr.fixed <- attributes(res.fixed)
int.fixed <-  attr.fixed$intercept
depvar <- as.character(attr.fixed$variables[2])
inddepvar.fixed <- attr.fixed$term.labels
inddepvar.fixed.nom <-inddepvar.fixed 
if(int.fixed > 0) inddepvar.fixed.nom <-c("intercept",inddepvar.fixed)

##########################################################
res.mixture <- terms(mixture)
attr.mixture <- attributes(res.mixture) 
int.mixture <-  attr.mixture$intercept  
inddepvar.mixture <- attr.mixture$term.labels 
inddepvar.mixture.nom <- inddepvar.mixture
if(int.mixture > 0) inddepvar.mixture.nom <- c("intercept",inddepvar.mixture)

##########################################################"
res.random <- terms(random)
attr.random <- attributes(res.random) 
int.random <-  attr.random$intercept  
inddepvar.random <- attr.random$term.labels 
inddepvar.random.nom <- inddepvar.random
if(int.random > 0) inddepvar.random.nom <- c("intercept",inddepvar.random)
##########################################################"
res.classmb <- terms(classmb)
attr.classmb <- attributes(res.classmb) 
int.classmb <-  attr.classmb$intercept  
inddepvar.classmb <- attr.classmb$term.labels 
inddepvar.classmb.nom <- inddepvar.classmb
inddepvar.classmb.nom <- c("intercept",inddepvar.classmb)


###########################################################
# intercept is always in inddepvar.classmb
var.exp <- unique(c(inddepvar.fixed,inddepvar.mixture,inddepvar.random,inddepvar.classmb))

nom.fixed <- inddepvar.fixed.nom
nom.mixture <- inddepvar.mixture.nom  

if(!(all(nom.mixture %in% nom.fixed))) stop("The covariates in mixture should be also included in the argument fixed")


Y0 <- data[,depvar] 
X0 <- as.data.frame(data[,var.exp])
names(X0) <- var.exp
if((any(is.na(X0))==TRUE)|(any(is.na(Y0))==TRUE))stop("The data should not contain any missing value")
 
n <- dim(data)[1]
if ((int.fixed+int.random)>0) X0<- cbind(intercept=rep(1,n),X0)
nom.X0 <- names(X0)
nvar.exp <- length(nom.X0)

IND <- data[,nom.subject]


#### INCLUSION PRIOR 
if(missing(prior)){ PRIOR <- seq(0,length=length(IND))} 
if(!missing(prior)){ 
PRIOR <- data[,nom.prior]
PRIOR[(is.na(PRIOR))] <- 0
}
####

ng0 <- ng
idiag0 <- as.integer(idiag)
nwg0 <- as.integer(nwg)

idea0 <- rep(0,nvar.exp)
idprob0 <- rep(0,nvar.exp)
idg0 <- rep(0,nvar.exp)

for (i in 1:nvar.exp)    {
 idea0[i] <- nom.X0[i]%in%inddepvar.random.nom
 idprob0[i] <- nom.X0[i]%in%inddepvar.classmb.nom      
 if(nom.X0[i]%in%nom.fixed & !(nom.X0[i]%in%nom.mixture)) idg0[i] <- 1 
 if(nom.X0[i]%in%nom.fixed & nom.X0[i]%in%nom.mixture) idg0[i] <- 2  
 }

if((int.fixed+int.random)>0) idprob0[1] <- 0



# on ordonne les donn es suivants la variable IND
matYX <- cbind(IND,PRIOR,Y0,X0)
matYXord <- matYX[sort.list(matYX[,1]),]
Y0 <- matYXord[,3]  
X0 <- matYXord[,-c(1,2,3)]
IND <- matYXord[,1]


#### INCLUSION PRIOR 
PRIOR <- matYXord[,2]
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

#-------------------------------------------------------------------------------
#definition du vecteur de parametre + initialisation
#-------------------------------------------------------------------------------
#####cas 1 : ng=1
b<-NULL
b1 <- NULL
NPROB <- 0
if(ng0==1| missing(B)){
NEF<-sum(idg0!=0)
b1[1:NEF]<-0
if(int.fixed > 0)  b1[1]<-mean(Y0)

if(idiag0==1){
NVC<-sum(idea0==1)
b1[(NEF+1):(NEF+NVC)]<-1}

if(idiag0==0){
kk<-sum(idea0==1) 
NVC<-(kk*(kk+1))/2
indice<-cumsum(1:kk)
bidiag<-rep(0,NVC)
bidiag[indice]<-1
b1[(NEF+1):(NEF+NVC)]<-bidiag
}
b1[NEF+NVC+1]<-1
NPM<-length(b1)
NW<-0
V <- rep(0,NPM*(NPM+1)/2) 
}

#####cas 2 : ng>=2
if(ng0>1){
NPROB<-(sum(idprob0==1)+1)*(ng0-1)
b[1:NPROB]<-0
NEF<-sum(idg0==1)+(sum(idg0==2))*ng0
if(idiag0==1)NVC<-sum(idea0==1)
if(idiag0==0){
kk<-sum(idea0==1) 
NVC<-(kk*(kk+1))/2}
NW<-nwg0*(ng0-1)
if(NW>0) b[(NPROB+NEF+NVC+1):(NPROB+NEF+NVC+NW)]<-1
NPM<-NPROB+NEF+NVC+NW+1
V <- rep(0,NPM*(NPM+3)/2)
} 



if(missing(B)){

if(ng0>1){
idea2 <- idea0
idprob2 <- rep(0,nv0)  
idg2 <- rep(0,nv0) 
idg2[idg0!=0] <- 1
NEF2<-sum(idg2==1)
NPM2<-NEF2+NVC+1
nwg2<-0
ng2<-1
ppi2<- rep(0,ns0)
pred_m_g2 <- rep(0,nobs0)
pred_ss_g2 <- rep(0,nobs0)

V2 <- rep(0,NPM2*(NPM2+1)/2)
best <- rep(0,NPM2)
init <- .Fortran("hetmixlin",as.double(Y0),as.double(X0),as.integer(prior2),as.integer(idprob2),as.integer(idea2),as.integer(idg2),as.integer(ns0),as.integer(ng2),as.integer(nv0),as.integer(nobs0),as.integer(nea0),as.integer(nmes0),as.integer(idiag0),as.integer(nwg2),npm=as.integer(NPM2),best=as.double(b1),V=as.double(V2),loglik=as.double(loglik),niter=as.integer(ni),conv=as.integer(istop),gconv=as.double(gconv),ppi2=as.double(ppi2),resid_m=as.double(resid_m),resid_ss=as.double(resid_ss),pred_m_g=as.double(pred_m_g2),pred_ss_g=as.double(pred_ss_g2),predRE=as.double(predRE),as.double(convB),as.double(convL),as.double(convG),as.integer(Maxiter),PACKAGE="lcmm")

k <- NPROB
l <- 0
t<- 0
for (i in 1:nvar.exp)    {
if(idg0[i]==1){
l <- l+1
t <- t+1
b[k+t] <- init$best[l]
}
if(idg0[i]==2){
	l <- l+1
	for (g in 1:ng){
	t <- t+1
	b[k+t] <- init$best[l]+(g-(ng+1)/2)*sqrt(init$V[l*(l+1)/2])
	}
}
}
b[(NPROB+NEF+1):(NPROB+NEF+NVC)] <-init$best[(NEF2+1):(NEF2+NVC)]
b[NPROB+NEF+NVC+NW+1] <- init$best[NPM2]
} 
if(ng0==1 ){
b <- b1
}
} else {if(length(B)!=NPM)stop("The length of the vector B is not correct")
 else {b <-B}
} 

se <- rep(0,length(b))
#------------------------------------------
#------nom au vecteur best
#--------------------------------------------

if(ng0==2)names(b)[1:NPROB]<-inddepvar.classmb.nom

if(ng0>2){
nom <-rep(inddepvar.classmb.nom,each=ng0-1)
nom1 <- paste(nom,c(1:(ng0-1)))
names(b)[1:NPROB]<-nom1
}

if(ng0==1) names(b)[1:NEF] <- inddepvar.fixed.nom

if(ng0>1){
nom1<- NULL
for (i in 1:nvar.exp) {
if(idg0[i]==2){ nom <- paste(nom.X0[i],c(1:ng0))
nom1 <- cbind(nom1,t(nom))}
if(idg0[i]==1) nom1 <- cbind(nom1,nom.X0[i])
}
names(b)[(NPROB+1):(NPROB+NEF)]<- nom1
}
if(NVC!=0)names(b)[(NPROB+NEF+1):(NPROB+NEF+NVC)] <- paste("varcov",c(1:(NVC)))
if(NW!=0)names(b)[(NPROB+NEF+NVC+1):(NPROB+NEF+NVC+NW)] <- paste("varprop",c(1:(ng0-1)))
names(b)[(NPROB+NEF+NVC+NW+1):(NPM)] <- "stderr"

N <- NULL
N[1] <- NPROB
N[2] <- NEF
N[3] <- NVC
N[4] <- NW

idiag <- as.integer(idiag0)
idea <- as.integer(idea0)
nv <- as.integer(nv0)


################ Sortie ###########################

out <- .Fortran("hetmixlin",as.double(Y0),as.double(X0),as.integer(prior0),as.integer(idprob0),as.integer(idea0),as.integer(idg0),as.integer(ns0),as.integer(ng0),as.integer(nv0),as.integer(nobs0),as.integer(nea0),as.integer(nmes0),as.integer(idiag0),as.integer(nwg0),as.integer(NPM),best=as.double(b),V=as.double(V),loglik=as.double(loglik),niter=as.integer(ni),conv=as.integer(istop),gconv=as.double(gconv),ppi2=as.double(ppi0),resid_m=as.double(resid_m),resid_ss=as.double(resid_ss),pred_m_g=as.double(pred_m_g),pred_ss_g=as.double(pred_ss_g),predRE=as.double(predRE),as.double(convB),as.double(convL),as.double(convG),as.integer(Maxiter),PACKAGE="lcmm")


### Creation du vecteur cholesky
Cholesky <- rep(0,(nea0*(nea0+1)/2))
if(idiag0==0){
Cholesky[1:NVC] <- out$best[(NPROB+NEF+1):(NPROB+NEF+NVC)]
### Construction de la matrice U 
U <- matrix(0,nrow=nea0,ncol=nea0)
U[upper.tri(U,diag=TRUE)] <- Cholesky[1:NVC]
z <- t(U) %*% U
out$best[(NPROB+NEF+1):(NPROB+NEF+NVC)] <- z[upper.tri(z,diag=TRUE)]
}
if(idiag0==1){
id <- 1:nea0
indice <- rep(id+id*(id-1)/2)
Cholesky[indice] <- out$best[(NPROB+NEF+1):(NPROB+NEF+nea0)]
out$best[(NPROB+NEF+1):(NPROB+NEF+NVC)] <- out$best[(NPROB+NEF+1):(NPROB+NEF+NVC)]**2 
} 

####################################################

predRE <- matrix(out$predRE,ncol=nea0,byrow=T)
predRE <- cbind(INDuniq,predRE)
colnames(predRE) <- c(nom.subject,inddepvar.random.nom)

ppi<- matrix(out$ppi2,ncol=ng0,byrow=TRUE)
classif<-apply(ppi,1,which.max)
ppi<-cbind(INDuniq,classif,ppi)
temp<-paste("prob",1:ng0,sep="")
colnames(ppi) <- c(nom.subject,"class",temp)
rownames(ppi) <- 1:ns0

pred_m_g <- matrix(out$pred_m_g,nrow=nobs0)
pred_ss_g <- matrix(out$pred_ss_g,nrow=nobs0)
pred_m <- Y0-out$resid_m
pred_ss <- Y0-out$resid_ss
pred <- cbind(IND,pred_m,out$resid_m,pred_ss,out$resid_ss,Y0,pred_m_g,pred_ss_g)

temp<-paste("pred_m",1:ng0,sep="")
temp1<-paste("pred_ss",1:ng0,sep="")
colnames(pred)<-c(nom.subject,"pred_m","resid_m","pred_ss","resid_ss","obs",temp,temp1) 

names(out$best)<-names(b)
btest <- out$best[1:length(inddepvar.fixed.nom)]
names(btest) <-inddepvar.fixed.nom
res <-list(ns=ns0,ng=ng0,idea0=idea0,idprob0=idprob0,idg0=idg0,loglik=out$loglik,best=out$best,V=out$V,gconv=out$gconv,conv=out$conv,call=cl,niter=out$niter,dataset=args$data,N=N,name.mat.cov=inddepvar.random.nom,idiag=idiag0,pred=pred,pprob=ppi,predRE=predRE,Xnames=nom.X0,cholesky=Cholesky)
class(res) <-"hlme"
#class(res) <-"lcmm"  
res
}
