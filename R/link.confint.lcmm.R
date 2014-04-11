

link.confint.lcmm <- function(x,ndraws=2000,Yvalues,plot=TRUE,...)
{
 ## verification des arguments
 if(missing(x)) stop("The model should be specified.")
 if(!(inherits(x,"lcmm"))) stop("To use only with \"lcmm\" objects")
 if(!missing(Yvalues) & x$linktype==3) warning("With thresholds links, no \"Yvalues\" should be specified. Default values will be used. \n")

 
 if(x$conv %in% c(1,2)) 
 {
  if(missing(Yvalues) | x$linktype==3)
  {
   new.transf <- FALSE
   Yvalues <- x$estimlink[,1]
  }
  else
  {
   new.transf <- TRUE
   Yvalues <- na.omit(Yvalues)
   if(any(Yvalues<x$estimlink[1,1]) | any(Yvalues>x$estimlink[nrow(x$estimlink),1])) stop("The values specified in \"Yvalues\" are not in the range of the outcome")
   Yvalues <- sort(Yvalues)
  }
  
  
  ##preparation des arguments pour Fortran
  npm <- length(x$best)
  best <- x$best
  if(x$idiag==0 & x$N[3]>0) best[(x$N[1]+x$N[2]+1):(x$N[1]+x$N[2]+x$N[3])] <- x$cholesky
  if(x$idiag==1 & x$N[3]>0) best[(x$N[1]+x$N[2]+1):(x$N[1]+x$N[2]+x$N[3])] <- sqrt(best[(x$N[1]+x$N[2]+1):(x$N[1]+x$N[2]+x$N[3])])
 
  if(x$linktype==0) ntrtot <- 2
  if(x$linktype==1) ntrtot <- 4
  if(x$linktype==2) ntrtot <- length(x$linknodes)+2
  if(x$linktype==3) ntrtot <- sum(x$ide)
  
  imoins <- x$N[1]+x$N[2]+x$N[3]+x$N[4]
  zitr <- x$linknodes
  maxnbzitr <- ifelse(x$linktype==2,length(x$linknodes),2)
  epsY <- x$epsY
  minY <- x$estimlink[1,1]
  maxY <- x$estimlink[nrow(x$estimlink),1]
  ny <- 1
  nsim <- length(Yvalues)
  
  if(x$linktype==3)
  {
   ide <- x$ide
   dimide <- length(ide)
  }
  else
  {
   ide <- rep(0,1) #pas utilise si pas ordinal
   dimide <- 1
  }
 
  ndraws <- as.integer(ndraws)
 
  Mat <- matrix(0,ncol=npm,nrow=npm)
   # que la partie sup utilisée donc OK si rien en bas
  Mat[upper.tri(Mat,diag=TRUE)]<- x$V
  Chol <- chol(Mat)
  Chol <- t(Chol)
 
 
 
  ## calcul des valeurs trasnformees si necessaire
  if(isTRUE(new.transf))
  {
   resFortran <- rep(0,nsim)
  
   out0 <- .Fortran("calculus_transfo",as.double(best),as.integer(npm),as.integer(ny),as.integer(x$linktype),
   as.integer(ntrtot),as.integer(imoins),as.double(zitr),as.integer(maxnbzitr),as.double(Yvalues),as.integer(nsim),
   as.double(minY),as.double(maxY),as.double(epsY),as.integer(ide),as.integer(dimide),transfo=as.double(resFortran),PACKAGE="lcmm")
   
   transfY <- out0$transfo
  }
  else
  {
   transfY <- x$estimlink[,2]
  }
  
  if(x$conv==1)
  {
   ## boucle pour b=1,...B :
   Hydraws <- NULL
   for (j in 1:ndraws)
   {
    bdraw <- rnorm(npm)
    bdraw <- best + Chol %*% bdraw
  
    resFortran <- rep(0,nsim)
    
    out <- .Fortran("calculus_transfo",as.double(bdraw),as.integer(npm),as.integer(ny),as.integer(x$linktype),
    as.integer(ntrtot),as.integer(imoins),as.double(zitr),as.integer(maxnbzitr),as.double(Yvalues),as.integer(nsim),
    as.double(minY),as.double(maxY),as.double(epsY),as.integer(ide),as.integer(dimide),transfo=as.double(resFortran),PACKAGE="lcmm")
    
    Hydraws <- cbind(Hydraws,out$transfo)
   }
   
   ## calcul des bornes IC
   f <- function(x)
   {
    quantile(x[!is.na(x)],probs=c(0.025,0.975))
   }
   
   Hydistr <- apply(Hydraws,1,FUN=f)
   borne_inf <- as.vector(Hydistr[1,])
   borne_sup <- as.vector(Hydistr[2,])
  
  
   ## resultat a renvoyer
   res <- data.frame(Yvalues=Yvalues,transfY=transfY,borne_inf=borne_inf,borne_sup=borne_sup)
  }
  
  if(x$conv==2)
  {
   borne_inf <- rep(NA,length(Yvalues))
   borne_sup <- rep(NA,length(Yvalues))
   res <- data.frame(Yvalues=Yvalues,transfY=transfY,borne_inf=borne_inf,borne_sup=borne_sup)
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
    if(x$linktype!=3)
    {    
     dots <- c(dots,list(xlim=range(c(transfY,borne_inf,borne_sup),na.rm=TRUE)))
    }
    else
    {
     dots <- c(dots,list(xlim=range(c(transfY[-c(1,length(transfY))],borne_inf[-c(1,length(transfY))],borne_sup[-c(1,length(transfY))]),na.rm=TRUE)))
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
    dots <- c(dots,list(col=1))
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
  
    copiex <- x
    ysim <- matrix(Yvalues,nsim,ny)
    transfo <- matrix(borne_sup,nsim,ny)
    estimlink <- as.vector(rbind(ysim,transfo))
    copiex$estimlink <- matrix(estimlink,nsim,2*ny)
    
    #tracer la borne sup
    do.call("plot.linkfunction",c(dots.bornes,list(x=copiex,add=TRUE,legend=NULL)))
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
