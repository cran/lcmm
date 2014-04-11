

link.confint.multlcmm <- function(x,ndraws=2000,Yvalues,plot=TRUE,...)
{
 ## verification des arguments
 if(missing(x)) stop("The model should be specified.")
 if(!(inherits(x,"multlcmm"))) stop("To use only with \"multlcmm\" objects")
 #if(!missing(Yvalues) & x$linktype==3) warning("With thresholds links, no \"Yvalues\" should be specified. Default values will be used. \n")


 if(x$conv %in% c(1,2))
 {  
  ny <- x$N[8]

  if(missing(Yvalues))
  {
   new.transf <- FALSE
   Yvalues <- as.vector(x$estimlink[,2*1:ny-1])
  }
  else
  {
   new.transf <- TRUE
   na.fail(Yvalues)
   Yvalues <- apply(Yvalues,2,sort)
   
   #controler si minY<Yvalues<maxY
   for(yk in 1:ny)
   {
    if(any(Yvalues[,yk]<x$estimlink[1,2*yk-1]) | any(Yvalues[,yk]>x$estimlink[nrow(x$estimlink),2*yk-1])) stop("The values specified in \"Yvalues\" are not in the range of the outcome")
   }
   
   Yvalues <- as.vector(Yvalues)
  }
 
 
  ##preparation des arguments pour Fortran
  npm <- length(x$best)
  best <- x$best
  if(x$N[4]>0)
  {
   if(x$idiag==0) best[x$N[3]+1:x$N[4]] <- x$cholesky[-1]
   else best[x$N[3]+1:x$N[4]] <- x$cholesky[c((1:(x$N[4]+1)*2:(x$N[4]+2))/2)[-1]]
  }
  
  ntrtot <- rep(NA,x$N[8])
  numSPL <- 0
  for (yk in 1:x$N[8])
  {
   if (x$linktype[yk]==0)
   {
    ntrtot[yk] <- 2
   }
   if (x$linktype[yk]==1)
   {
   ntrtot[yk] <- 4
   }
   if (x$linktype[yk]==2)
   {
    numSPL <-  numSPL+1
    ntrtot[yk] <- x$nbnodes[numSPL]+2
   }
  }
 
  imoins <- sum(x$N[3:8])
  zitr <- x$linknodes
  nbzitr <- rep(2,ny)
  nbzitr[which(x$linktype==2)] <- x$nbnodes
  maxnbzitr <- max(nbzitr)
  epsY <- x$epsY
  minY <- x$estimlink[1,2*1:ny-1]
  maxY <- x$estimlink[nrow(x$estimlink),2*1:ny-1]
 
  nsim <- length(Yvalues)/ny
 
  ide <- matrix(0,nrow=1,ncol=ny) #pas encore de threshold link
  dimide <- rep(1,ny)
 
  ndraws <- as.integer(ndraws)
 
  Mat <- matrix(0,ncol=npm,nrow=npm)
   # que la partie sup utilisée donc OK si rien en bas
  Mat[upper.tri(Mat,diag=TRUE)]<- x$V
  Chol <- chol(Mat)
  Chol <- t(Chol)
 
 
 
  ## calcul des valeurs trasnformees si necessaire
  if(isTRUE(new.transf))
  {
   resFortran <- rep(0,nsim*ny)
 
   out0 <- .Fortran("calculus_transfo",as.double(best),as.integer(npm),as.integer(ny),as.integer(x$linktype),
   as.integer(ntrtot),as.integer(imoins),as.double(zitr),as.integer(maxnbzitr),as.double(Yvalues),as.integer(nsim),
   as.double(minY),as.double(maxY),as.double(epsY),as.integer(ide),as.integer(dimide),transfo=as.double(resFortran),PACKAGE="lcmm")
 
   transfY <- out0$transfo
  }
  else
  {
   transfY <- as.vector(x$estimlink[,2*1:ny])
  }
 
  if(x$conv==1)
  {
   ## boucle pour b=1,...B :
   Hydraws <- NULL
   for (j in 1:ndraws)
   {
    bdraw <- rnorm(npm)
    bdraw <- best + Chol %*% bdraw
  
    resFortran <- rep(0,nsim*ny)
  
    out <- .Fortran("calculus_transfo",as.double(bdraw),as.integer(npm),as.integer(ny),as.integer(x$linktype),
    as.integer(ntrtot),as.integer(imoins),as.double(zitr),as.integer(maxnbzitr),as.double(Yvalues),as.integer(nsim),
    as.double(minY),as.double(maxY),as.double(epsY),as.integer(ide),as.integer(dimide),transfo=as.double(resFortran),PACKAGE="lcmm")
  
    Hydraws <- cbind(Hydraws,out$transfo)
   }
  
   ## calcul des bornes IC
   f <- function(x)
   {
    quantile(x[!is.na(x)],probs=c(0.025,0.5,0.975))
   }
  
   Hydistr <- apply(Hydraws,1,FUN=f)
   transfY <- as.vector(Hydistr[2,])
   borne_inf <- as.vector(Hydistr[1,])
   borne_sup <- as.vector(Hydistr[3,])
   
   ## resultat a renvoyer
   Yname <- rep(x$Ynames,each=nsim)
   res <- data.frame(Yname=Yname,Yvalues=Yvalues,transfY=transfY,borne_inf=borne_inf,borne_sup=borne_sup)    
  }
  
  if(x$conv==2)
  {
   ## resultat a renvoyer
   Yname <- rep(x$Ynames,each=nsim)
   borne_inf <- rep(NA,length(Yname))
   borne_sup <- rep(NA,length(Yname))
   res <- data.frame(Yname=Yname,Yvalues=Yvalues,transfY=transfY,borne_inf=borne_inf,borne_sup=borne_sup)   
  }
  
   
  ## graphique
  if(isTRUE(plot))
  {
   copiex <- x
   ysim <- matrix(Yvalues,nsim,ny)
   transfo <- matrix(transfY,nsim,ny)
   estimlink <- as.vector(rbind(ysim,transfo))
   copiex$estimlink <- matrix(estimlink,nsim,2*ny)
   
   dots <- list(...)

   if(length(list(...)$xlim)==0)
   {
    if(all(x$linktype!=3))
    {    
     dots <- c(dots,list(xlim=range(c(transfY,borne_inf,borne_sup),na.rm=TRUE)))
    }
    else
    {
     big <- c(nsim*(which(x$linktype==3)-1)+1,nsim*which(x$linktype==3))
     dots <- c(dots,list(xlim=range(c(transfY[-big],borne_inf[-big],borne_sup[-big]),na.rm=TRUE)))
    }
   }   

      
   #tracer la transfo
   do.call("plot.linkfunction",c(dots,list(x=copiex)))
 

   if(length(list(...)$lwd)==0)
   {
    dots <- c(dots,list(lwd=1))
   }

   if(length(list(...)$col)==0)
   {
    dots <- c(dots,list(col=rainbow(ny)))
   }
    
 
   if(x$conv==1)
   {
    #on ne garde que les options lwd et col et on impose lty=2
    dots.bornes <- dots[intersect(names(dots),c("lwd","col"))]
    dots.bornes <- c(dots.bornes,lty=2) 
   
    copiex <- x
    ysim <- matrix(Yvalues,nsim,ny)
    transfo <- matrix(borne_inf,nsim,ny)
    estimlink <- as.vector(rbind(ysim,transfo))
    copiex$estimlink <- matrix(estimlink,nsim,2*ny)
    
    #tracer la borne inf
    do.call("plot.linkfunction",c(dots.bornes,list(x=copiex,add=TRUE,legend=NULL)))
    #lcmm::plot.linkfunction(copiex,add=TRUE,legend=NULL,lty=2,lwd=lwd1,col=col1)
  
    copiex <- x
    ysim <- matrix(Yvalues,nsim,ny)
    transfo <- matrix(borne_sup,nsim,ny)
    estimlink <- as.vector(rbind(ysim,transfo))
    copiex$estimlink <- matrix(estimlink,nsim,2*ny)
    
    #tracer la borne sup
    do.call("plot.linkfunction",c(dots.bornes,list(x=copiex,add=TRUE,legend=NULL)))
    #lcmm::plot.linkfunction(copiex,add=TRUE,legend=NULL,lty=2,lwd=lwd1,col=col1)
   }
  }
 }
 else
 {
  cat("Output can not be produced since the program stopped abnormally.")
  res <- NA 
 }

 return(res)
}

link.confint <- function(x,ndraws=2000,Yvalues,plot=TRUE,...) UseMethod("link.confint")
