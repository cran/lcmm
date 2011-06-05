
############################### Fonction Jointlcmm ###################################

Jointlcmm <-
function(fixed,mixture,random,subject,classmb,ng=1,idiag=FALSE,nwg=FALSE,survival,hazard="Weibull",hazardtype="Specific",hazardnodes=NULL,TimeDepVar=NULL,data,B,convB=0.0001,convL=0.0001,convG=0.0001,maxiter=500,nsim=100,prior){


cl <- match.call()
m <- match.call(expand = FALSE)
    m$fixed <- m$mixture <- m$random <- m$subject <- m$classmb <-m$ng<-m$idiag<-m$nwg<-m$B<-m$convB <- m$convL <-m$convG<-m$prior<-m$maxiter<-m$hazardnodes<-m$hazard<-m$hazardtype<-m$TimeDepVar<-m$nsim<-m$... <-NULL
args <- as.list(match.call(Jointlcmm))[-1]

#nom.subject <- as.character(args$subject)
nom.subject <- as.character(subject)

#### INCLUSION PRIOR
nom.prior <- as.character(args$prior)

#### ERROR MESSAGES
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
if(class(survival)!="formula") stop("The argument survival must be a formula")

if(missing(data)){ stop("The argument data should be specified and defined as a data.frame")} 
if(missing(subject)){ stop("The argument subject must be specified in any model even without random-effects")} 


###########################################################
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
    
    m<-model.frame(formula(m$survival),data=data)
    
    depvar1 <- model.extract(m, "response")

    if (!inherits(depvar1, "Surv")){ 
        stop("Response must be a survival object")
    }
### for covariates in survival
### nom des variables a droite dans survie (avec et sans mixture)
    inddep.surv <- attr(res.evt, "term.labels")

### ind.mixture = toutes les var avec mixture dans survie
    ind.mixture <- untangle.specials(res.evt, "mixture", 1)

### variables avec mixture
    dcomp <- ind.mixture$vars
    dcomp <- gsub("mixture\\(","",dcomp)
    dcomp <- gsub("\\)","",dcomp)
    inddepvar.Mixt  <- dcomp


 
### inddepvar.survival = toutes les variables avec pas de mixture dans survie
    inddepvar.noMixt <- inddep.surv[!(inddep.surv %in% ind.mixture$vars)]

### ind.mixture$vars et ajouter a inddepvar.survival
    inddepvar.survival <- c(inddepvar.noMixt,inddepvar.Mixt)

    mt <- attr(m, "terms")  
    X <- if (!is.empty.model(mt)){ 
        model.matrix(mt, m, contrasts)
    }
################# type du modele de survie
    attr.surv.type <- attr(depvar1,"type")

    if (attr.surv.type=="right"){
        idtrunc0 <-0
        Tevent <- depvar1[,1]
        Devent <- depvar1[,2]
    }
    if (attr.surv.type=="right"){
        Tentry <- rep(0,length(Tevent))
    }

    if (attr.surv.type=="counting"){
        idtrunc0 <- 1
        Tentry <- depvar1[,1]
        Tevent <- depvar1[,2]
        Devent <- depvar1[,3]
    }

    if (max(Devent)==1) evt <- 1
    if (max(Devent)==0) evt <- 0

if (!(attr.surv.type %in% c("right","counting"))){
   stop("Jointlcmm handles only 'right' and 'counting' types of Survival data")
}


##############   COVARIATES       ##########################
# intercept is always in inddepvar.classmb
var.exp <- unique(c(inddepvar.fixed,inddepvar.mixture,inddepvar.random,inddepvar.classmb,inddepvar.survival))

nom.fixed <- inddepvar.fixed.nom
nom.mixture <- inddepvar.mixture.nom  
if(!(all(nom.mixture %in% nom.fixed))) stop("The covariates in mixture should be also included in the argument fixed")



##############   DATA      ##########################
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



#### DEFINITION INDICATORS 
ng0 <- ng
idiag0 <- as.integer(idiag)
nwg0 <- as.integer(nwg)
idea0 <- rep(0,nvar.exp)
idprob0 <- rep(0,nvar.exp)
idg0 <- rep(0,nvar.exp)
idxevt <- rep(0,nvar.exp)

for (i in 1:nvar.exp){
 idea0[i] <- nom.X0[i]%in%inddepvar.random.nom
 idprob0[i] <- nom.X0[i]%in%inddepvar.classmb.nom      
 if(nom.X0[i]%in%nom.fixed & !(nom.X0[i]%in%nom.mixture)) idg0[i] <- 1 
 if(nom.X0[i]%in%nom.fixed & nom.X0[i]%in%nom.mixture) idg0[i] <- 2 
 if((nom.X0[i]%in%inddepvar.survival) & !(nom.X0[i]%in%inddepvar.Mixt)) idxevt[i] <- 1
 if((nom.X0[i]%in%inddepvar.survival) & (nom.X0[i]%in%inddepvar.Mixt)) idxevt[i] <- 2 
 }
 
if((int.fixed+int.random)>0) idprob0[1] <- 0



###### DATA SORTING on IND variable
matYX <- cbind(IND,PRIOR,Y0,X0)
matYXord <- matYX[sort.list(matYX[,1]),]
Y0 <- matYXord[,3]
X0 <- matYXord[,-c(1,2,3)]
IND <- matYXord[,1]
PRIOR <- matYXord[,2]
PRIOR <-as.integer(as.vector(PRIOR))

X0<-as.numeric(as.matrix(X0))
Y0<-as.numeric(as.matrix(Y0))
nmes0<-as.vector(table(IND))

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

#### Parametre timeDepVar  ####

nvdepsurv <- 0
tsurvint <- tsurv
indsurvint <- rep(0,ns0)

if(all.equal(TimeDepVar,NULL)==F){

   nvdepsurv <- 1
   tsurvint <- tsurv
   TimeDepVar[(is.na(TimeDepVar))] <- max(tsurv)
   indsurvint[TimeDepVar<tsurv] <- 1
   tsurvint[TimeDepVar<tsurv] <- TimeDepVar[TimeDepVar<tsurv]
}




##### hazard specification

# revoir le test ci-dessous : peut-être inclu dans un autre
     carac <- strsplit(hazard,split="-")
     carac <- unlist(carac)   
     if(!(length(carac) %in% c(1,2,3))){stop("Please check and revise the hazard argument according to the format specified in the help.")}
     
# test sur hazardtype
     if (!(hazardtype %in% c("PH","Common","Specific"))) stop("Only 'Specific', 'PH' or 'Common' hazardtype can be specified corresponding respectively to class-specific hazards, hazards proportional over latent classes or hazards common over classes")
     
     haz <-strsplit(hazard,split="")
     haz <- unlist(haz)
     haz <- grep("-",haz)

# revoir ce test car pourrait y avoir que splines ou piecewise	
     if((all.equal(length(haz),0)==T)==T){
	if(!(hazard %in% c("Weibull","piecewise","splines"))){
		stop("Only 'Weibull', 'piecewise' or 'splines' hazard can be specified in hazard argument")
	} 
     }else{
### contrÃ´le du separateur de l'argument hazard

	if(!(length(haz) %in% c(1,2))) stop("With splines or piecewise baseline function, the separator of hazard argument must be only '-'")  
	
	test <- strsplit(hazard,split="-")
	test <- unlist(test)   
	
	if(all.equal(length(test),2)==T){
		if(!all(test[1:2] %in% c("splines","piecewise","equi","manual","quant"))){
			stop ("With splines or Piecewise baseline function, hazard argument should contain only 'splines','piecewise','equi','manual','quant'")
		}		
	}else{

		if(!all(test[2:3] %in% c("splines","piecewise","equi","manual","quant"))){
			stop ("With splines or piecewise baseline function, hazard argument should contain only 'splines','piecewise','equi','manual','quant'")
		}
	}	 
     }     
 

nrisqtot <- 0

### Pour Weibull
if(all.equal(hazard,"Weibull")==T){
       nz0=2       
       typrisq0 <- 2
       nprisq0 <- 2 
}
else{
## test si on entre les noeuds manuellement	

	if(("manual" %in% unlist(strsplit(hazard,split="-")))){	
		if(all.equal(length(unlist(strsplit(hazard,split="-"))),2)==T){
			nz0 <- length(hazardnodes)+2
			typrisq0 <- switch(test[3],"piecewise"=1,"splines"=3)
			nprisq0 <- switch(test[3],"piecewise"=nz0-1,"splines"=nz0+2)
		}else{
			nz0 <- as.integer(test[1])
			nz1 <- length(hazardnodes)+2
			if(!(all.equal(nz0,nz1)==T)){
#				cat("Warning: on garde nz1","\n") 
				nz0 <- nz1				
			}
			
			typrisq0 <- switch(test[3],"piecewise"=1,"splines"=3)
			nprisq0 <- switch(test[3],"piecewise"=nz0-1,"splines"=nz0+2)	
		}
	
	}else{
		nz0 <- as.integer(test[1])
		typrisq0 <- switch(test[3],"piecewise"=1,"splines"=3)
		nprisq0 <- switch(test[3],"piecewise"=nz0-1,"splines"=nz0+2)	
	
	}
}
risqcom0 <- switch(hazardtype,"Specific"=0,"PH"=2,"Common"=1)
if (evt != 0){   
  nrisqtot <-  switch(hazardtype,"Specific"=nprisq0*ng0,"PH"=nprisq0+ng0-1,"Common"=nprisq0)   
} 


##### parametre nvarxevt  	    
nvarxevt <- sum(idxevt==1)+(sum(idxevt==2))*ng0
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
zi0[1] <- min(c(tsurv0,tsurv))
zi0[nz0] <- max(c(tsurv0,tsurv))
}
if(idtrunc0==0){
zi0[1] <- min(c(tsurv))
zi0[nz0] <- max(c(tsurv))
}
minT <- zi0[1]
maxT <- zi0[nz0]
if((maxT-minT) < 0) stop("Please check the time of event variable. It seems that all the times are equal.")

if (!(length(grep("Weibull",hazard)) > 0)){
##### si "equi"
    if(nz0-2 <= 0){
         stop("Splines or piecewise baseline function should include at least 2 nodes (and at least 5 are recommended)")
    }else{
        if((all.equal("equi",test[2])==T)==T){
            pas=as.double(maxT-minT)/as.double(nz0-1)
            for(i in 2:(nz0-1)){
                zi0[i] <- zi0[i-1]+ pas
            }
        }
	
##### si "manual" 
if(all.equal(length(unlist(strsplit(hazard,split="-"))),2)==T){
        if((all.equal("manual",test[1])==T)==T){
            if (is.null(hazardnodes)){
                 stop("If 'manual' option is specified for the splines or piecewise baseline hazard function, hazardnodes argument should include the list of interior nodes")
            }else{
                zi0[2:(nz0-1)] <- hazardnodes[1:(nz0-2)]
            }
         }
}else{
        if((all.equal("manual",test[2])==T)==T){
            if (is.null(hazardnodes)){
                 stop("If 'manual' option is specified for the splines or piecewise baseline hazard function, hazardnodes argument should include the list of interior nodes")
            }else{
		    hazardnodes <- sort(hazardnodes)
		    zi0[2:(nz0-1)] <- hazardnodes[1:(nz0-2)]
            }
         }

}
##### si "quant"
       if((all.equal("quant",test[2])==T)==T){
          pas <-c(1:(nz0-2))/(nz0-1) 
          zi0 [2:(nz0-1)] <- quantile(TSURV,probs=pas)
       }
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
   NVARXEVTinit <- sum(idxevt!=0)

   b1[1:(nprisq0)]<-(1/nprisq0)
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
   b1[(NEF+NVC+1)]<-1
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
      b[(NEF+NVC+NW+1)]<-1
      NPM<-NEF+NVC+NW+1
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
         NPM2<-NEF2+NVC+1
         nwg2<-0
         ng2<-1
         ppi2 <- rep(0,ns0)
	 ppitest2 <- rep(0,ns0)
         pred_m_g2 <- rep(0,nobs0)
         pred_ss_g2 <- rep(0,nobs0)
         V2 <- rep(0,NPM2*(NPM2+1)/2)

#      cat("Binit",b1,"\n")

init <- .Fortran("Jointhet",as.double(Y0),as.double(X0),as.integer(initial),as.integer(idprob2),as.integer(idea2),as.integer(idg2),as.integer(idxevt2),as.integer(ns0),as.integer(ng2),as.integer(nv0),as.integer(nobs0),as.integer(nea0),as.integer(nmes0),as.integer(idiag0),as.integer(nwg2),as.integer(NPM2),best=as.double(b1),V=as.double(V2),as.double(loglik),as.integer(ni),as.integer(istop),as.double(gconv),as.double(ppi2),as.double(ppitest2),as.double(resid_m),as.double(resid_ss),as.double(pred_m_g2),as.double(pred_ss_g2),as.double(predRE),as.double(convB),as.double(convL),as.double(convG),as.integer(maxiter),as.integer(evt),as.integer(typrisq0),as.integer(idtrunc0),as.integer(risqcom0),as.integer(nz0),as.double(zi0),as.integer(nvdepsurv),as.double(tsurv0),as.double(tsurv),as.double(tsurvint),as.integer(devt),as.integer(indsurvint),as.double(statsc),as.double(risq_est),as.double(surv_est),as.integer(nsim),as.double(time_est))


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
	           b[(NPROB+nprisq0*(g-1)+1):(NPROB+nprisq0*g)] <- abs(init$best[1:nprisq0])+(g-(ng+1)/2)*sqrt(init$V[ident])
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
	         b[k+t] <- init$best[l]+(g-(ng+1)/2)*sqrt(init$V[l*(l+1)/2])
              }
           }
        }

## k pour le nouveau vecteur B (ie nrisqtot+nvarxevt)
## l pour le vecteur B issu de init (ie nprisq0+NVARXEVTinit)

        k <- NPROB+nrisqtot+nvarxevt
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
	         b[k+t] <- init$best[l]+(g-(ng+1)/2)*sqrt(init$V[l*(l+1)/2])
              }
           }
        }
        b[(NEF+1):(NEF+NVC)] <-init$best[(NEF2+1):(NEF2+NVC)]
        b[(NEF+NVC+NW+1)] <- init$best[NPM2]
     }
 
     if(ng0==1){
        b <- b1
     }
} else {
    if(length(B)!=NPM){
       cat("taille",length(B),"NPM",NPM)
       stop("The length of the vector B is not correct")
    } else{
    b <-B
    }
} 

#------------------------------------------
#------nom au vecteur best
#--------------------------------------------

if(ng0==2){
   names(b)[1:NPROB]<-inddepvar.classmb.nom
}



if(ng0>2){
   nom <-rep(inddepvar.classmb.nom,each=ng0-1)
   nom1 <-paste(nom,".",c(1:(ng0-1)),sep="")  
   names(b)[1:NPROB]<-nom1
}


if(ng0==1) {
    names(b)[(1):(nrisqtot)] <- paste("Survparm",c(1:(nprisq0)),sep="")
    names(b)[(nrisqtot+1):(nprisq0+nvarxevt)] <- inddepvar.survival#evt
    names(b)[(nrisqtot+NVARXEVTinit+1):(NEF)] <- inddepvar.fixed.nom 
}


if(ng0>1){
# survie
   if (risqcom0==1){
       names(b)[(NPROB+1):(NPROB+nrisqtot)] <- paste("Survparm",c(1:(nrisqtot)),sep="")
   }


   if (risqcom0==2) {
      names(b)[(NPROB+1):(NPROB+nprisq0)] <- paste("Survparm",c(1:(nprisq0)),sep="")
      names(b)[(NPROB+nprisq0+1):(NPROB+nrisqtot)] <- paste("SurvPH",c(1:(ng0-1)),sep="")
   }
   if (risqcom0==0){
      for (g in 1:ng0){
	names(b)[((NPROB+nprisq0*(g-1))+1):(NPROB+nprisq0*g)] <- paste("Survparm",c(1:(nprisq0)),".",g,sep="")
      }
   }
#noms evol
   nom1<- NULL
   nom2<- NULL
   for (i in 1:nvar.exp) {
      if(idg0[i]==2){ 
          nom <- paste(nom.X0[i],".",c(1:ng0),sep="")	  
          nom1 <- cbind(nom1,t(nom))
      }
      if(idg0[i]==1){
         nom1 <- cbind(nom1,nom.X0[i])
      }
# noms survie  
     
#idxevt0
      if(idxevt[i]==2){ 
          nom <- paste(nom.X0[i],".",c(1:ng0),sep="")
          nom2 <- cbind(nom2,t(nom))
      }
      if(idxevt[i]==1) nom2 <- cbind(nom2,nom.X0[i])
   }  
    names(b)[(nrisqtot+NPROB+1):(nrisqtot+nvarxevt+NPROB)]<- nom2
    names(b)[(nrisqtot+nvarxevt+NPROB+1):NEF]<- nom1

}


if(NVC!=0)names(b)[(NEF+1):(NEF+NVC)] <- paste("Varcov",c(1:(NVC)),sep="")
if(NW!=0)names(b)[(NEF+NVC+1):(NEF+NVC+NW)] <- paste("Varprop",c(1:(ng0-1)),sep="")
names(b)[(NEF+NVC+NW+1):(NPM)] <- "Stderr"


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
#############################
idiag <- as.integer(idiag0)
idea <- as.integer(idea0)
nv <- as.integer(nv0)


################ Sortie ###########################


################## Verification des parametres d'entree du modele ##################


#      cat("Binit",b,"\n")

cat("Be patient, Jointlcmm program is running ... \n")


out <- .Fortran("Jointhet",as.double(Y0),as.double(X0),as.integer(prior0),as.integer(idprob0),as.integer(idea0),as.integer(idg0),as.integer(idxevt),as.integer(ns0),as.integer(ng0),as.integer(nv0),as.integer(nobs0),as.integer(nea0),as.integer(nmes0),as.integer(idiag0),as.integer(nwg0),as.integer(NPM),best=as.double(b),V=as.double(V),loglik=as.double(loglik),niter=as.integer(ni),conv=as.integer(istop),gconv=as.double(gconv),ppi2=as.double(ppi0),ppitest2=as.double(ppitest0),resid_m=as.double(resid_m),resid_ss=as.double(resid_ss),pred_m_g=as.double(pred_m_g),pred_ss_g=as.double(pred_ss_g),predRE=as.double(predRE),as.double(convB),as.double(convL),as.double(convG),as.integer(maxiter),as.integer(evt),as.integer(typrisq0),as.integer(idtrunc0),as.integer(risqcom0),as.integer(nz0),as.double(zi0),as.integer(nvdepsurv),as.double(tsurv0),as.double(tsurv),as.double(tsurvint),as.integer(devt),as.integer(indsurvint),statsc=as.double(statsc),risq_est=as.double(risq_est),surv_est=as.double(surv_est),nsim=as.integer(nsim),time_est=as.double(time_est))


if (!(out$conv %in% c(1,2))){
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


### Creation du vecteur cholesky



Cholesky <- rep(NA,(nea0*(nea0+1)/2))
if(idiag0==0){
   Cholesky[1:NVC] <- out$best[(NEF+1):(NEF+NVC)]
### Construction de la matrice U 
   U <- matrix(0,nrow=nea0,ncol=nea0)
   U[upper.tri(U,diag=TRUE)] <- Cholesky[1:NVC]
   z <- t(U) %*% U
   out$best[(NEF+1):(NEF+NVC)] <- z[upper.tri(z,diag=TRUE)]
}
if(idiag0==1){
   id <- 1:nea0
   indice <- rep(id+id*(id-1)/2)
   Cholesky[indice] <- out$best[(NEF+1):(NEF+nea0)]
   out$best[(NEF+1):(NEF+NVC)] <- out$best[(NEF+1):(NEF+NVC)]**2 
} 


####################################################

predRE <- matrix(out$predRE,ncol=nea0,byrow=T)
predRE <- cbind(INDuniq,predRE)
colnames(predRE) <- c(nom.subject,inddepvar.random.nom)

ppi<- matrix(out$ppi2,ncol=ng0,byrow=TRUE)

if(!(out$conv %in% c(1,2))){
	classif <-rep(NA,ns0)
}else{
	classif<-apply(ppi,1,which.max)
}


ppi<-cbind(INDuniq,classif,ppi)
temp<-paste("probYT",1:ng0,sep="")
colnames(ppi) <- c(nom.subject,"class",temp)
rownames(ppi) <- 1:ns0

ppitest<- matrix(out$ppitest2,ncol=ng0,byrow=TRUE)
if(!(out$conv %in% c(1,2))){
	classif <-rep(NA,ns0)
}else{
	classif<-apply(ppitest,1,which.max)
}


ppitest<-cbind(INDuniq,classif,ppitest)
temp<-paste("probY",1:ng0,sep="")
colnames(ppitest) <- c(nom.subject,"class",temp)
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


pred <- cbind(IND,pred_m,out$resid_m,pred_ss,out$resid_ss,Y0,pred_m_g,pred_ss_g)


temp<-paste("pred_m",1:ng0,sep="")
temp1<-paste("pred_ss",1:ng0,sep="")
colnames(pred)<-c(nom.subject,"pred_m","resid_m","pred_ss","resid_ss","obs",temp,temp1) 

names(out$best)<-names(b)




surv_est <- matrix(out$surv_est,nrow=nsim)
risq_est <- matrix(out$risq_est,nrow=nsim)
predSurv <- cbind(out$time_est,risq_est,surv_est)

temp<-paste("RiskFct",1:ng0,sep="")
temp1<-paste("CumRiskFct",1:ng0,sep="")
colnames(predSurv)<-c("time",temp,temp1) 


# rajouter des arguments en sortie
############################# AJOUT survnodes le 10/03/2011
res <-list(ns=ns0,ng=ng0,idea0=idea0,idprob0=idprob0,idg0=idg0,loglik=out$loglik,best=out$best,V=out$V,gconv=out$gconv,conv=out$conv,call=cl,niter=out$niter,dataset=args$data,N=N,name.mat.cov=inddepvar.random.nom,idiag=idiag0,pred=pred,pprob=ppi,pprobY=ppitest,predRE=predRE,Xnames=nom.X0,cholesky=Cholesky,CIstat=out$statsc,predSurv=predSurv,hazard=list(typrisq0,hazardtype,zi0))

###########################
class(res) <-c("Jointlcmm")
res
}
