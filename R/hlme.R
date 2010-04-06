hlme <-
function(fixed,mixture,random,subject,classmb,ng=1,idiag=FALSE,nwg=FALSE,data,B){
cat("Be patient. The program is running ... \n")

cl <- match.call()
args <- as.list(match.call(hlme))[-1]

nom.subject <- as.character(args$subject)

if(!missing(mixture) & ng==1) stop("No mixture possible with ng=1")
if(missing(mixture) & ng>1) stop("The argument mixture should be specified for ng > 1")
if(!missing(classmb) & ng==1) stop("No classmb possible with ng=1")
if(missing(random)) random <- ~-1
if(missing(fixed)) stop("fixed must be specified by a formula")
if(missing(classmb)) classmb <- ~-1
if(missing(mixture)) mixture <- ~-1
if(ng==1&nwg==TRUE) stop("The argument nwg should be FALSE")


if(class(fixed)!="formula") stop("fixed must be a formula")
if(class(mixture)!="formula") stop("mixture must be a formula")
if(class(random)!="formula") stop("random must be a formula")
if(class(classmb)!="formula") stop("classmb must be a formula")

if(missing(data)){ stop("the data should be specified and must be a data.frame")} 
if(missing(subject)){ stop("the argument subject must be provided")} 



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
#if(int.classmb > 0) inddepvar.classmb.nom <- c("intercept",inddepvar.classmb)
inddepvar.classmb.nom <- c("intercept",inddepvar.classmb)


###########################################################
var.exp <- unique(c(inddepvar.fixed,inddepvar.mixture,inddepvar.random,inddepvar.classmb))

nom.fixed <- inddepvar.fixed.nom
nom.mixture <- inddepvar.mixture.nom  

if(!(all(nom.mixture %in% nom.fixed))) stop("the arguments mixture should be included in the argument fixed")


Y0 <- data[,depvar] 
X0 <- as.data.frame(data[,var.exp])
names(X0) <- var.exp
if((any(is.na(X0))==TRUE)|(any(is.na(Y0))==TRUE))stop("The data should not contained missing value")
 
n <- dim(data)[1]
if ((int.fixed+int.random)>0) X0<- cbind(intercept=rep(1,n),X0)
nom.X0 <- names(X0)
nvar.exp <- length(nom.X0)

IND <- data[,nom.subject]
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

# on ordonne les données suivants la variable IND
matYX <- cbind(IND,Y0,X0)
matYXord <- matYX[sort.list(matYX[,1]),]
Y0 <- matYXord[,2]  
X0 <- matYXord[,-c(1,2)]


X0<-as.numeric(as.matrix(X0))
Y0<-as.numeric(as.matrix(Y0))
nmes0<-as.vector(table(IND))
ns0<-length(nmes0)


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
NPM<-NPROB+NEF+NVC+NW+1}



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
se2<-rep(0,NPM2)
init <- .Fortran("hetmixlin",as.double(Y0),as.double(X0),as.integer(idprob2),as.integer(idea2),as.integer(idg2),as.integer(ns0),as.integer(ng2),as.integer(nv0),as.integer(nobs0),as.integer(nmes0),as.integer(idiag0),as.integer(nwg2),npm=as.integer(NPM2),best=as.double(b1),se=as.double(se2),loglik=as.double(loglik),niter=as.integer(ni),conv=as.integer(istop),gconv=as.double(gconv),ppi2=as.double(ppi2),resid_m=as.double(resid_m),resid_ss=as.double(resid_ss),pred_m_g=as.double(pred_m_g2),pred_ss_g=as.double(pred_ss_g2),PACKAGE="lcmm")
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
	b[k+t] <- init$best[l]+(g-(ng+1)/2)*init$se[l]
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

N <- NULL
N[1] <- NPROB
N[2] <- NEF
N[3] <- NVC
N[4] <- NW


################## CHANGEMENTS TAILLE ####################"


out <- .Fortran("hetmixlin",as.double(Y0),as.double(X0),as.integer(idprob0),as.integer(idea0),as.integer(idg0),as.integer(ns0),as.integer(ng0),as.integer(nv0),as.integer(nobs0),as.integer(nmes0),as.integer(idiag0),as.integer(nwg0),npm=as.integer(NPM),best=as.double(b),se=as.double(se),loglik=as.double(loglik),niter=as.integer(ni),conv=as.integer(istop),gconv=as.double(gconv),ppi2=as.double(ppi0),resid_m=as.double(resid_m),resid_ss=as.double(resid_ss),pred_m_g=as.double(pred_m_g),pred_ss_g=as.double(pred_ss_g),PACKAGE="lcmm")


ppi<- matrix(out$ppi2,ncol=ng0,byrow=TRUE)
classif<-apply(ppi,1,which.max)
ppi<-cbind(classif,ppi)
temp<-paste("prob",1:ng0,sep="")
colnames(ppi) <- c("class",temp)
rownames(ppi) <- 1:ns0

pred_m_g <- matrix(out$pred_m_g,nrow=nobs0)
pred_ss_g <- matrix(out$pred_ss_g,nrow=nobs0)
pred_m <- Y0-out$resid_m
pred_ss <- Y0-out$resid_ss
pred <- cbind(pred_m,out$resid_m,pred_ss,out$resid_ss,Y0,pred_m_g,pred_ss_g)

temp<-paste("pred_m",1:ng0,sep="")
temp1<-paste("pred_ss",1:ng0,sep="")
colnames(pred)<-c("pred_m","resid_m","pred_ss","resid_ss","obs",temp,temp1) 

names(out$best)<-names(b)
btest <- out$best[1:length(inddepvar.fixed.nom)]
names(btest) <-inddepvar.fixed.nom
res <-list(ns=ns0,ng=ng0,idea0=idea0,idprob0=idprob0,idg0=idg0,loglik=out$loglik,best=out$best,se=out$se,gconv=out$gconv,conv=out$conv,call=cl,niter=out$niter,dataset=args$data,N=N,name.mat.cov=inddepvar.random.nom,idiag=idiag0,pred=pred,pprob=ppi,Xnames=nom.X0)
class(res) <-"hlme"  
res

}

