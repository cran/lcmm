
multlcmm <- function(fixed,mixture,random,subject,classmb,ng=1,idiag=FALSE,nwg=FALSE,randomY=FALSE,link="linear",intnodes=NULL,epsY=0.5,cor=NULL,data,B,convB=0.0001,convL=0.0001,convG=0.0001,maxiter=100,verbose=FALSE,nsim=100,prior,range=NULL,na.action=1)
{
 ptm<-proc.time()
 cat("Be patient, multlcmm is running ... \n")

 cl <- match.call()

 nom.subject <- as.character(subject)

 #### INCLUSION PRIOR
 nom.prior <- NULL
 if(!missing(prior)) nom.prior <- as.character(prior)
 ####

 if(!missing(mixture) & ng==1) stop("No mixture can be specified with ng=1")
 if(missing(mixture) & ng>1) stop("The argument mixture has to be specified for ng > 1")
 if(!missing(classmb) & ng==1) stop("No classmb can be specified with ng=1")
 if(missing(random)) stop("At least one random effect is required")
 if(random==~-1) stop("At least one random effect is required")
 if(missing(fixed)) stop("The argument Fixed must be specified in any model")
 if(missing(classmb) & ng==1) classmb <- ~-1
 if(missing(classmb) & ng>1) classmb <- ~1
 if(missing(mixture)) mixture <- ~-1
 if(ng==1&nwg==TRUE) stop ("The argument nwg should be FALSE for ng=1")


 if(class(fixed)!="formula") stop("The argument fixed must be a formula")
 if(class(mixture)!="formula") stop("The argument mixture must be a formula")
 if(class(random)!="formula") stop("The argument random must be a formula")
 if(class(classmb)!="formula") stop("The argument classmb must be a formula")
 if(missing(data)){ stop("The argument data should be specified and defined as a data.frame")}
 if(nrow(data)==0) stop("Data should not be empty")
 if(missing(subject)){ stop("The argument subject must be specified in any model even without random-effects")}
 if(any(link=="thresholds"))  stop("The link function thresholds is not available in multivariate case")
 if(all(link %in% c("linear","beta")) & !is.null(intnodes)) stop("Intnodes should only be specified with splines links")

 if(!(na.action%in%c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument")

 ### test de l'argument cor
 ncor0 <- 0
 cor.type <- cl$cor[1]
 cor.time <- cl$cor[2]
 cor <- paste(cor.type,cor.time,sep="-")
 if (!isTRUE(all.equal(cor,character(0))))
 {
  if (substr(cor,1,2)=="AR") { ncor0 <- 2 }
  else if (substr(cor,1,2)=="BM") { ncor0 <- 1  }
  else { stop("The argument cor must be of type AR or BM") }

  if(!(strsplit(cor,"-")[[1]][2] %in% colnames(data))) stop("Unaible to find time variable from argument cor in data")
  else { cor.var.time <- strsplit(cor,"-")[[1]][2] }
 }
 ### fin test argument cor

 #pour acces aux attributs des formules
 afixed <- terms(fixed, specials=c("factor","contrast"))
 if(attr(afixed,"intercept")==0) stop("An intercept should appear in fixed for identifiability purposes")
 amixture <- terms(mixture, specials=c("factor"))
 if (any(!(all.vars(amixture) %in% all.vars(afixed)))) stop("Variables in mixture should also appear in fixed")
 arandom <- terms(random, specials=c("factor"))
 aclassmb <- terms(classmb, specials=c("factor"))
  #fixed sans contrast
 fixed2 <- gsub("contrast","",fixed)
 fixed2 <- formula(paste(fixed2[2],fixed2[3],sep="~"))
 afixed2 <- terms(fixed2)
  #contrast
 contr <- ~-1
 if(!is.null(attr(afixed,"specials")$contrast))
 {
  vcontr <- attr(afixed,"term.labels")[setdiff(attr(afixed,"specials")$contrast-1,untangle.specials(afixed,"contrast",2)$terms)]
  vcontr <- gsub("contrast","",vcontr)
  contr <- as.formula(paste("~-1+",paste(vcontr,collapse="+")))
 }
 acontr <- terms(contr)
 
 
 #tjrs intercept dans classmb
 if(attr(aclassmb,"intercept")==0 & ng>1)
 {
  attr(aclassmb,"intercept") <- 1
  cat("The formula in classmb should always include an intercept. An intercept has been added.")
 }

 ###liste des outcomes
 nomsY <- as.character(attr(afixed,"variables")[2])
 nomsY <- strsplit(nomsY,split=" + ",fixed=TRUE)
 nomsY <- as.vector(nomsY[[1]])
 ny0 <- length(nomsY)

 #pas de contrast ni randomY si un seul Y
 if(ny0<2 & length(attr(afixed,"specials")$contrast)) stop("No contrast can be included with less than two outcomes")
 if(ny0<2 & randomY==TRUE) stop("With less than 2 outcomes randomY should be FALSE")

 ###liste des variables utilisees  (sans les interactions et sans les Y)
 ttesLesVar <- colnames(get_all_vars(afixed,data=data[1,]))
 ttesLesVar <- c(ttesLesVar, colnames(get_all_vars(amixture,data=data[1,])))
 ttesLesVar <- c(ttesLesVar, colnames(get_all_vars(arandom,data=data[1,])))
 ttesLesVar <- c(ttesLesVar, colnames(get_all_vars(aclassmb,data=data[1,])))
 if (ncor0>0) ttesLesVar <- unique(c(ttesLesVar,cor.var.time))
 else ttesLesVar <- unique(ttesLesVar)
 ttesLesVar <- setdiff(ttesLesVar, nomsY)

 ###subset de data avec les variables utilisees
 newdata <- data[,c(nom.subject,nomsY,ttesLesVar,nom.prior)]
 if(!is.null(nom.prior))
 {
  prior <- newdata[,nom.prior]
  newdata[which(is.na(prior)),nom.prior] <- 0
 }

 ###un data frame par outcome et creation Y0
 dataY <- paste("data",nomsY,sep="")
 Y0 <- NULL
 IND <- NULL
 outcome <- NULL
 prior <- NULL
 data0 <- NULL
 for (k in 1:ny0)
 {
  dtemp <- newdata[,c(nom.subject,nomsY[k],ttesLesVar,nom.prior)]
  #enlever les NA
  linesNA <- apply(dtemp,2,function(v) which(is.na(v)))
  linesNA <- unique(unlist(linesNA))
  if(na.action==1 & length(linesNA)>0) dtemp <- dtemp[-linesNA,]
  if(na.action==2 & length(linesNA)>0) stop("Data contains missing values")
  assign(dataY[k],dtemp)
  Y0 <- c(Y0, dtemp[,nomsY[k]])
  IND <- c(IND, dtemp[,nom.subject])
  outcome <- c(outcome,rep(nomsY[k],nrow(dtemp)))
  if(!is.null(nom.prior)) prior <- c(prior, dtemp[,nom.prior])
  data0 <- rbind(data0, dtemp[,setdiff(colnames(dtemp),nomsY[k])])   #dataset sans NA avec les covariables utilisees; obs ordonnees par outcome
 }

 ###prior=0 si pas specifie
 if(is.null(prior)) prior <- rep(0,length(Y0))

 ###creation de X0 (ttes les var + interactions)

 Xfixed <- model.matrix(fixed2[-2], data=data0)
 Xmixture <- model.matrix(mixture, data=data0)
 Xrandom <- model.matrix(random, data=data0)
 Xclassmb <- model.matrix(classmb, data=data0)
 Xcontr <- model.matrix(contr,data=data0)

 #changer les noms des interactions
 tousLesTermes <- unique(c(labels(afixed2),labels(arandom),labels(amixture),labels(aclassmb),labels(acontr)))
 MatVarTerm <- apply(matrix(tousLesTermes,ncol=1),1,function(x) ttesLesVar %in% unlist(strsplit(x,":",fixed=TRUE)) + 0)
 rownames(MatVarTerm) <- ttesLesVar
 colnames(MatVarTerm) <- tousLesTermes

 uMat <- unique(MatVarTerm,MARGIN=2)
 dbleMat <- MatVarTerm[,setdiff(colnames(MatVarTerm),colnames(uMat)),drop=FALSE]

 if(ncol(dbleMat))
 {
  for(i in 1:ncol(dbleMat))
  {
   icol <- which(apply(uMat,2,all.equal,current=dbleMat[,i])=="TRUE")
   aGarder <- colnames(uMat[,icol,drop=FALSE])
   aRempl <- colnames(dbleMat[,i,drop=FALSE])

   colnames(Xfixed)[which(colnames(Xfixed)==aRempl)] <- aGarder
   colnames(Xmixture)[which(colnames(Xmixture)==aRempl)] <- aGarder
   colnames(Xrandom)[which(colnames(Xrandom)==aRempl)] <- aGarder
   colnames(Xclassmb)[which(colnames(Xclassmb)==aRempl)] <- aGarder
   colnames(Xcontr)[which(colnames(Xcontr)==aRempl)] <- aGarder
  }
 }



 X0 <- cbind(Xfixed, Xrandom, Xclassmb, Xmixture, Xcontr)
 nom.unique <- unique(colnames(X0))
 X0 <- X0[,nom.unique]

 if (ncor0>0)
 {
  if(!(cor.var.time %in% colnames(X0)))
  {
   X0 <- cbind(X0, data0[,cor.var.time])
   colnames(X0) <- c(nom.unique, cor.var.time)
  }
 }

 #X0 <- as.data.frame(X0)
 X0 <- as.matrix(X0)
 ###X0 fini


 ###test de link
 if (length(link)!=1 & length(link)!=ny0) stop("One link per outcome should be specified")
 if(length(link)==1 & ny0>1)
 {
  link <- rep(link, ny0)
  intnodes <- rep(intnodes, ny0)
 }

 idlink0 <- rep(2,ny0)
 idlink0[which(link=="linear")] <- 0
 idlink0[which(link=="beta")] <- 1

 if (any(idlink0==1))
 {
  if (epsY<=0)
  {
   epsY <- 0.5
   cat("Argument 'epsY' should be a definite positive real. It is changed to the default value of 0.5. \n")
  }
 }

 spl <- strsplit(link[which(idlink0==2)],"-")
 if(any(sapply(spl,length)!=3)) stop("Invalid argument link")

 nySPL <- length(spl)
 nybeta <- sum(idlink0==1)
 #remplir range si pas specifie
 if(!is.null(range) & length(range)!=2*(nySPL+nybeta)) stop("Length of vector range is not correct.")
 if((length(range)==2*(nySPL+nybeta)) & (nySPL+nybeta>0))
 {
  ind12 <- which(idlink0==1 | idlink0==2)
  for (i in 1:(nySPL+nybeta))
  {
   rg <- range(get(dataY[ind12[i]]))
   if(rg[1]>range[2*(i-1)+1] | rg[2]<range[2*(i-1)+2]) stop("The range specified do not cover the entire range of the data")
  }
 }
 if((is.null(range) & (nybeta+nySPL)>0) | length(range)!=2*(nySPL+nybeta))
 {
  range <- NULL
  for(k in which(idlink0!=0))
  {
   range <- c(range, min(get(dataY[k])[,nomsY[k]]), max(get(dataY[k])[,nomsY[k]]))
  }
 }



 nbzitr0 <- rep(2,ny0) #nbzitr0 = nb de noeuds si splines, 2 sinon
 nbnodes <- NULL  #que pour les splines
 spltype <- NULL
 if(nySPL>0)
 {
 for (i in 1:nySPL)
 {
  nbnodes <- c(nbnodes, spl[[i]][1])
  spltype <- c(spltype, spl[[i]][2])
  if(spl[[i]][3] != "splines") stop("Invalid argument link")
 }
 }
 nbnodes <- as.numeric(nbnodes)
 nbzitr0[which(idlink0==2)] <- nbnodes

 #test splines
 if(!(all(spltype %in% c("equi","quant","manual")))) stop("The location of the nodes should be 'equi', 'quant' or 'manual'")

 #intnodes2 : contient tous les noeuds interieurs (pas seulement ceux de manual)
 intnodes2 <- rep(NA,sum(nbnodes-2))
 nb <- 0
 nbspl <- 0
 for (k in 1:ny0)
 {
  if (idlink0[k]!=2) next
  else
  {
   nbspl <- nbspl+1

   if(spltype[nbspl]=="manual")
   {
    nodes <- intnodes[(nb+1):(nb+nbnodes[nbspl]-2)]
    intnodes2[(sum(nbnodes[1:nbspl]-2)-(nbnodes[nbspl]-2)+1):sum(nbnodes[1:nbspl]-2)] <-  nodes
    nb <- nb+nbnodes[nbspl]-2

    if(any(nodes <= range[2*(nbspl-1)+1]) | any(nodes >= range[2*nbspl])) stop("Interior nodes must be in the range of the outcome")
   }

    if(spltype[nbspl]=="equi")
    {
     nodes <- seq(range[2*(nbspl-1)+1], range[2*nbspl], length.out=nbnodes[nbspl])
     nodes <- nodes[-nbnodes[nbspl]]
     nodes <- nodes[-1]
     intnodes2[(sum(nbnodes[1:nbspl]-2)-(nbnodes[nbspl]-2)+1):sum(nbnodes[1:nbspl]-2)] <- nodes
    }

    if(spltype[nbspl]=="quant")
    {
     nodes <- quantile(get(dataY[k])[,nomsY[k]], probs=seq(0,1,length.out=nbnodes[nbspl]))
     if(length(unique(nodes)) != length(nodes)) stop(paste("Some nodes are equal for link number",k,"; Please try to reduce the number of nodes or use manual location."))
     nodes <- nodes[-nbnodes[nbspl]]
     nodes <- nodes[-1]
     intnodes2[(sum(nbnodes[1:nbspl]-2)-(nbnodes[nbspl]-2)+1):sum(nbnodes[1:nbspl]-2)] <- as.vector(nodes)
    }

  }
 }

 if(nb != length(intnodes)) stop("The length of intnodes is not correct")

 #remplir zitr
 m <- 0
 if(nySPL>0) m <- max(nbnodes)
 zitr <- matrix(0,max(m,2),ny0)
 nb12 <- 0
 nbspl <- 0
 for (k in 1:ny0)
 {
  if(idlink0[k]==0) zitr[1:2,k] <- c(min(get(dataY[k])[,nomsY[k]]),max(get(dataY[k])[,nomsY[k]]))

  if(idlink0[k]==1)
  {
   nb12 <- nb12 + 1
   zitr[1:2,k] <- range[2*(nb12-1)+1:2]
  }

  if(idlink0[k]==2)
  {
   nb12 <- nb12+1
   nbspl <- nbspl+1
   zitr[2:(nbzitr0[k]-1),k] <- intnodes2[ifelse(nbspl==1,0,1)*sum(nbnodes[1:(nbspl-1)]-2) + 1:(nbnodes[nbspl]-2)]
   zitr[1,k] <- range[2*(nb12-1)+1]
   zitr[nbnodes[nbspl],k]  <- range[2*nb12]
  }
 }

 ###uniqueY0 et indiceY0
 uniqueY0 <- NULL
 indiceY0 <- NULL
 nvalSPL0 <- NULL
 nb <- 0
 for (k in 1:ny0)
 {
  if(idlink0[k]!=2)
  {
   indiceY0 <- c(indiceY0, rep(0,length(get(dataY[k])[,nomsY[k]])))
   next
  }

  yk <- get(dataY[k])[,nomsY[k]]
  uniqueTemp <- sort(unique(yk))
  permut <- order(order(yk))  # sort(y)[order(order(y))] = y
  indice <- rep(1:length(uniqueTemp), as.vector(table(yk)))
  indiceTemp <- nb + indice[permut]

  nb <- nb + length(uniqueTemp)
  uniqueY0 <- c(uniqueY0, uniqueTemp)
  indiceY0 <- c(indiceY0, indiceTemp)
  nvalSPL0 <- c(nvalSPL0, length(uniqueTemp))
 }


 ###ordonner les mesures par individu
 IDnum <- as.numeric(IND)
 matYX <- cbind(IDnum,IND,prior,Y0,indiceY0,outcome,X0)
 matYXord <- matYX[order(IDnum),]
 Y0 <- matYXord[,4]
 X0 <- matYXord[,-c(1,2,3,4,5,6)]
 #X0 <- as.matrix(X0)  a remettre si X0 <- as.data.frame(X0) remis l.157
 IDnum <- matYXord[,1]
 IND <- matYXord[,2]
 outcome <- matYXord[,6]
 indiceY0 <- matYXord[,5]
 prior0 <- matYXord[,3]


 ###parametres pour hetmixContMult
 ns0 <- length(unique(IND))
 ng0 <- ng
 nv0 <- dim(X0)[2]
 nobs0 <- length(Y0)
 idiag0 <- ifelse(idiag==TRUE,1,0)
 nwg0 <- ifelse(nwg==TRUE,1,0)
 nalea0 <- ifelse(randomY==TRUE,ny0,0)

 loglik <- 0
 ni <- 0
 istop <- 0
 gconv <- rep(0,3)
 ppi0 <- rep(0,ns0*ng0)
 resid_m <- rep(0,nobs0)
 resid_ss <- rep(0,nobs0)
 pred_m_g <- rep(0,nobs0*ng0)
 pred_ss_g <- rep(0,nobs0*ng0)
 Yobs <- rep(0,nobs0)
 #predRE <- rep(0,ns0*nea0)
 predRE_Y <- rep(0,ns0*nalea0)
 rlindiv <- rep(0,ns0)
 marker <- rep(0,nsim*ny0)
 transfY <- rep(0,nsim*ny0)
 Ydiscrete <- 0
 UACV <- 0
 vraisdiscret <- 0


 ###nmes0
 nmes0 <- matrix(0,ns0,ny0)
 for (k in 1:ny0)
 {
  INDpresents <- which(unique(IND) %in% get(dataY[k])[,nom.subject])
  nmes0[INDpresents,k] <- as.vector(table(get(dataY[k])[,nom.subject]))
 }


 ###remplir idprob, etc
 idprob0 <- colnames(X0) %in% colnames(Xclassmb) +0

 idea0 <- colnames(X0) %in% colnames(Xrandom) +0

 idg0 <- (colnames(X0) %in% colnames(Xfixed)) + (colnames(X0) %in% colnames(Xmixture))+0

 if (ncor0>0)idcor0 <- colnames(X0) %in% cor.var.time +0
 else idcor0 <- rep(0,nv0)

 idcontr0 <- colnames(X0) %in% colnames(Xcontr) +0

 if(any(idcontr0==1 & idg0==0)) stop("Variables in contrast should also appear in fixed")

 nea0 <- sum(idea0)
 predRE <- rep(0,ns0*nea0)

 #nombre total de parametres
 NPM <- (ng0-1)*sum(idprob0) + sum(idg0==1)-1 + ng0*sum(idg0==2) + ncor0 + (ny0-1)*sum(idcontr0) +
        ifelse(idiag0==1,nea0,nea0*(nea0+1)/2)-1 + (ng0-1)*nwg0 + nalea0 + ny0 +
        2*sum(idlink0==0) + 4*sum(idlink0==1) + sum(nbnodes+2)

 V <- rep(0, NPM*(NPM+1)/2)  #pr variance des parametres

 nef <- (ng0-1)*sum(idprob0) + sum(idg0==1)-1 + ng0*sum(idg0==2) + (ny0-1)*sum(idcontr0)
 ncontr <- (ny0-1)*sum(idcontr0)
 nvc <- ifelse(idiag0==1,nea0,nea0*(nea0+1)/2)-1
 nw <- (ng0-1)*nwg0
 ntrtot0 <- nbzitr0+2
 ntrtot0[which(idlink0==0)] <- 2
 nprob <- sum(idprob0)*(ng0-1)

 #valeurs initiales
 if(!(missing(B)))
 {
  if (length(B)==NPM) b <- B
  else stop(paste("Vector B should be of length",NPM))
 }
 else
 {
  b <- rep(0,NPM)
  if (nvc>0)
  {
   if(idiag==1) b[nef+1:nvc] <- rep(1,nvc)
   if(idiag==0)
   {
    init.nvc <- diag(nea0)
    init.nvc <- init.nvc[upper.tri(init.nvc, diag=TRUE)]
    b[nef+1:nvc] <- init.nvc[-1]
   }
  }
  if(nwg0>0) b[nef+nvc+1:nw] <- 1
  if(ncor0==1) b[nef+nvc+nw+1] <- 1
  if(ncor0==2) b[nef+nvc+nw+1:2] <- c(0,1)

  b[nef+nvc+nw+ncor0+1:ny0] <-  1

  if(nalea0>0) b[nef+nvc+nw+ncor0+ny0+1:nalea0] <- 1

  for(k in 1:ny0)
  {
   if(idlink0[k]==0)
   {
    b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])-1] <- mean(get(dataY[k])[,nomsY[k]])
    b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])] <- 1
   }
   if(idlink0[k]==1)
   {
    b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])-3] <- 0
    b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])-2] <- -log(2)
    b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])-1] <- 0.7
    b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])] <- 0.1
   }
   if(idlink0[k]==2)
   {
    b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])-ntrtot0[k]+1] <- -2
    b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])-ntrtot0[k]+2:ntrtot0[k]] <- 0.1
   }
  }
 }

 #------------------------------------------
 #------nom au vecteur best
 #--------------------------------------------

 nom.X0 <- colnames(X0)
 nom.X0[nom.X0=="(Intercept)"] <- "intercept"
 if(ng0>=2)
 {
  nom <-rep(nom.X0[idprob0==1],each=ng0-1)
  nom1 <- paste(nom," class",c(1:(ng0-1)),sep="")
  names(b)[1:nprob]<-nom1
 }


 if(ng0==1) names(b)[1:(nef-ncontr)] <- nom.X0[-1][idg0[-1]!=0]
 if(ng0>1){
 	nom1<- NULL
 	for (i in 1:nv0) {
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
 names(b)[(nprob+1):(nef-ncontr)]<- nom1
 }

 if(idlink0[1]==0) names(b)[nef+nvc+nw+ncor0+ny0+nalea0+1:ntrtot0[1]]<- c("Linear 1","Linear 2")
 if(idlink0[1]==1) names(b)[nef+nvc+nw+ncor0+ny0+nalea0+1:ntrtot0[1]]<- paste("Beta",c(1:ntrtot0[1]),sep="")
 if(idlink0[1]==2) names(b)[nef+nvc+nw+ncor0+ny0+nalea0+1:ntrtot0[1]]<- paste("I-splines",c(1:ntrtot0[1]),sep="")
 if(ny0>1)
 {
 for (yk in 2:ny0)
 {
  if(idlink0[yk]==0) names(b)[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:(yk-1)])+1:ntrtot0[yk]]<- c("Linear 1","Linear 2")
  if(idlink0[yk]==1) names(b)[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:(yk-1)])+1:ntrtot0[yk]]<- paste("Beta",c(1:ntrtot0[yk]),sep="")
  if(idlink0[yk]==2) names(b)[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:(yk-1)])+1:ntrtot0[yk]]<- paste("I-splines",c(1:ntrtot0[yk]),sep="")
 }
 }
 if(nvc!=0)names(b)[nef+1:nvc] <- paste("varcov",c(1:nvc))
 if(nw!=0)names(b)[nef+nvc+1:nw] <- paste("varprop class",c(1:(ng0-1)))

 names(b)[nef+nvc+nw+ncor0+1:ny0] <- paste("std.err",1:ny0)

 if(ncor0>0) {names(b)[nef+nvc+nw+1:ncor0] <- paste("cor",1:ncor0,sep="")}
 if(nalea0!=0) names(b)[nef+nvc+nw+ncor0+ny0+1:nalea0] <- paste("std.randomY",1:ny0,sep="")
 if(ncontr!=0) names(b)[(nef-ncontr+1):nef] <- paste("contrast",paste(rep(1:sum(idcontr0),each=ny0-1),rep(1:(ny0-1),sum(idcontr0)),sep=""),sep="")


 #initialisation pour ng>1
 if(missing(B) & ng0>1)
 {
  prior02 <- rep(0,nobs0)
  idprob02 <- rep(0,nv0)
  idg02 <- idg0
  idg02[idg02==2] <- 1
  idcontr02 <- rep(0,nv0)
  ng02 <- 1
  nwg02 <- 0
  nalea02 <- 0

  nef2 <- sum(idg0!=0)-1
  NPM2 <- nef2+ nvc+ncor0+ny0+sum(ntrtot0)
  b2 <- c(n=rep(0,nef2),b[(nef+1):NPM])
  ind1 <- which(substr(names(b2),1,7)=="varprop")
  if(length(ind1)) b2 <- b2[-ind1]
  ind2 <- which(substr(names(b2),1,11)=="std.randomY")
  if(length(ind2)) b2 <- b2[-ind2]

  V2 <- rep(0,NPM2*(NPM2+1)/2)
  loglik2 <- 0
  ppi02 <- rep(0,ns0)
  pred_m_g2 <- rep(0,nobs0)
  pred_ss_g2 <- rep(0,nobs0)
  predRE_Y2 <- rep(0,ns0)
  maxiter2 <- min(75,maxiter)
  convB2 <- max(0.01,convB)
  convL2 <- max(0.01,convL)
  convG2 <- max(0.01,convG)
  verbose2 <- FALSE

  if (verbose) cat("Computing initial values... \n")
  init <- .Fortran("hetmixContMult",as.double(Y0),as.double(X0),as.integer(prior02),as.integer(idprob02),as.integer(idea0),
  as.integer(idg02),as.integer(idcor0),as.integer(idcontr02),as.integer(ny0),as.integer(ns0),as.integer(ng02),
  as.integer(nv0),as.integer(nobs0),as.integer(nea0),as.integer(nmes0),as.integer(idiag0),as.integer(nwg02),
  as.integer(ncor0),as.integer(nalea02),as.integer(NPM2),best=as.double(b2),V=as.double(V2),loglik=as.double(loglik2),
  niter=as.integer(ni),conv=as.integer(istop),gconv=as.double(gconv),ppi2=as.double(ppi02),resid_m=as.double(resid_m),
  resid_ss=as.double(resid_ss),pred_m_g=as.double(pred_m_g2),pred_ss_g=as.double(pred_ss_g2),predRE=as.double(predRE),
  predRE_Y=as.double(predRE_Y2),as.double(convB2),as.double(convL2),as.double(convG2),as.integer(maxiter2),as.integer(verbose2),as.double(epsY),
  as.integer(idlink0),as.integer(nbzitr0),as.double(zitr),as.double(uniqueY0),as.integer(indiceY0),as.integer(nvalSPL0),
  marker=as.double(marker),transfY=as.double(transfY),as.integer(nsim),Yobs=as.double(Yobs),as.integer(Ydiscrete),
  vraisdiscret=as.double(vraisdiscret),UACV=as.double(UACV),rlindiv=as.double(rlindiv),PACKAGE="lcmm")


  l <- 0
  t <- 0
  for (i in 1:nv0)
  {
   if(idg0[i]==1 & i>1)
   {
    l <- l+1
    t <- t+1
    b[nprob+t] <- init$best[l]
   }
   if(idg0[i]==2)
   {
    if (i==1)
    {
	   for (g in 2:ng0)
     {
      t <- t+1
	    b[nprob+t] <- - 0.5*(g-1)
     }
    }
    if (i>1)
    {
	   l <- l+1
	   for (g in 1:ng0)
     {
	    t <- t+1
	    if(init$conv==1) b[nprob+t] <- init$best[l]+(g-(ng0+1)/2)*sqrt(init$V[l*(l+1)/2])
	    else b[nprob+t] <- init$best[l]+(g-(ng0+1)/2)*init$best[l]
	   }
    }
   }
  }

  if(nvc>0) b[nef+1:nvc] <-init$best[nef2+1:nvc]
  if (ncor0>0) {b[nef+nvc+nw+1:ncor0] <- init$best[nef2+nvc+1:ncor0]}
  b[nef+nvc+nw+ncor0+1:ny0] <- init$best[nef2+nvc+ncor0+1:ny0]
  b[(nef+nvc+nw+ncor0+ny0+nalea0+1):NPM] <-init$best[(nef2+nvc+ncor0+ny0+1):NPM2]

  if(verbose==TRUE) cat("initial parameters : \n",b,"\n")
 }


 #estimation
 out <- .Fortran("hetmixContMult",as.double(Y0),as.double(X0),as.integer(prior0),as.integer(idprob0),as.integer(idea0),
 as.integer(idg0),as.integer(idcor0),as.integer(idcontr0),as.integer(ny0),as.integer(ns0),as.integer(ng0),
 as.integer(nv0),as.integer(nobs0),as.integer(nea0),as.integer(nmes0),as.integer(idiag0),as.integer(nwg0),
 as.integer(ncor0),as.integer(nalea0),as.integer(NPM),best=as.double(b),V=as.double(V),loglik=as.double(loglik),
 niter=as.integer(ni),conv=as.integer(istop),gconv=as.double(gconv),ppi2=as.double(ppi0),resid_m=as.double(resid_m),
 resid_ss=as.double(resid_ss),pred_m_g=as.double(pred_m_g),pred_ss_g=as.double(pred_ss_g),predRE=as.double(predRE),
 predRE_Y=as.double(predRE_Y),as.double(convB),as.double(convL),as.double(convG),as.integer(maxiter),as.integer(verbose),as.double(epsY),
 as.integer(idlink0),as.integer(nbzitr0),as.double(zitr),as.double(uniqueY0),as.integer(indiceY0),as.integer(nvalSPL0),
 marker=as.double(marker),transfY=as.double(transfY),as.integer(nsim),Yobs=as.double(Yobs),as.integer(Ydiscrete),
 vraisdiscret=as.double(vraisdiscret),UACV=as.double(UACV),rlindiv=as.double(rlindiv),PACKAGE="lcmm")


 ### Creation du vecteur cholesky
 Cholesky <- rep(0,(nea0*(nea0+1)/2))
 if(idiag0==0 & nvc>0){
 Cholesky[1:(nvc+1)] <- c(1,out$best[nef+1:nvc])
 # Construction de la matrice U
 U <- matrix(0,nrow=nea0,ncol=nea0)
 U[upper.tri(U,diag=TRUE)] <- Cholesky[1:(nvc+1)]
 z <- t(U) %*% U
 out$best[nef+1:nvc] <- z[upper.tri(z,diag=TRUE)][-1]
 }
 if(idiag0==1 & nvc>0){
 id <- 1:nea0
 indice <- rep(id+id*(id-1)/2)
 Cholesky[indice] <- c(1,out$best[nef+1:nvc])
 out$best[nef+1:nvc] <- out$best[nef+1:nvc]**2
 }

 ###predictions
 predRE <- matrix(out$predRE,ncol=nea0,byrow=T)
 predRE <- data.frame(unique(IND),predRE)
 colnames(predRE) <- c(nom.subject,nom.X0[idea0!=0])

 if (nalea0!=0)
 {
  predRE_Y <- matrix(out$predRE_Y,ncol=ny0,byrow=TRUE)
  predRE_Y <- data.frame(unique(IND),predRE_Y)
  colnames(predRE_Y)  <- c(nom.subject,nomsY)
 }
 else
 {
  predRE_Y <- rep(NA,nalea0*ns0)
 }

 ###ppi
 if(ng0>1) {
 ppi<- matrix(out$ppi2,ncol=ng0,byrow=TRUE)
 }
 else {
 ppi <- matrix(rep(1,ns0),ncol=ng0)
 }

 classif<-apply(ppi,1,which.max)
 ppi<-data.frame(unique(IND),classif,ppi)
 temp<-paste("prob",1:ng0,sep="")
 colnames(ppi) <- c(nom.subject,"class",temp)
 rownames(ppi) <- 1:ns0


 ###pred
 pred_m_g <- matrix(out$pred_m_g,nrow=nobs0)
 pred_ss_g <- matrix(out$pred_ss_g,nrow=nobs0)
 pred_m <- out$Yobs-out$resid_m
 pred_ss <- out$Yobs - out$resid_ss
 pred <- data.frame(IND,outcome,pred_m,out$resid_m,pred_ss,out$resid_ss,out$Yobs,pred_m_g,pred_ss_g)

 temp<-paste("pred_m",1:ng0,sep="")
 temp1<-paste("pred_ss",1:ng0,sep="")
 colnames(pred)<-c(nom.subject,"Yname","pred_m","resid_m","pred_ss","resid_ss","obs",temp,temp1)
 rownames(pred) <- NULL

 ###estimlink
 ysim <- matrix(out$marker,nsim,ny0)
 transfo <- matrix(out$transfY,nsim,ny0)
 estimlink <- as.vector(rbind(ysim,transfo))
 estimlink <- matrix(estimlink,nsim,2*ny0)
 colnames(estimlink) <- paste(c("","transf"),rep(nomsY, each=2),sep="")


 N <- NULL
 N[1] <- (ng0-1)*sum(idprob0)
 N[2] <- (ny0-1)*sum(idcontr0)
 N[3] <- (ng0-1)*sum(idprob0) + sum(idg0==1)-1 + ng0*sum(idg0==2) + (ny0-1)*sum(idcontr0)  #nef
 N[4] <- ifelse(idiag0==1,nea0,nea0*(nea0+1)/2)-1  #nvc
 N[5] <- (ng0-1)*nwg0
 N[6] <- nalea0
 N[7] <- ncor0
 N[8] <- ny0
 N[9] <- nobs0

 nom.X0[nom.X0=="(Intercept)"] <- "Intercept"

 res <-list(ns=ns0,ng=ng0,idea0=idea0,idprob0=idprob0,idg0=idg0,idcontr0=idcontr0,idcor0=idcor0,loglik=out$loglik,best=out$best,V=out$V,gconv=out$gconv,conv=out$conv,call=cl,niter=out$niter,N=N,idiag=idiag0,pred=pred,pprob=ppi,predRE=predRE,predRE_Y=predRE_Y,Ynames=nomsY,Xnames=nom.X0,Xnames2=ttesLesVar,cholesky=Cholesky,estimlink=estimlink,epsY=epsY,linktype=idlink0,linknodes=zitr,nbnodes=nbnodes)
 names(res$best) <- names(b)
 class(res) <-c("multlcmm")

 cost<-proc.time()-ptm
 cat("The program took", round(cost[3]), "seconds \n")

 res
}
