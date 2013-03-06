predictY.multlcmm <- function(x,newdata,methInteg=0,nsim=20,draws=FALSE,ndraws=2000,na.action=1,...)
{
if(missing(newdata)) stop("The argument newdata should be specified")
if(missing(x)) stop("The argument x should be specified")
if (!inherits(x, "multlcmm")) stop("use only with \"lcmm\" or \"multlcmm\" objects")
if (!all(x$Xnames2 %in% colnames(newdata))) stop(paste(c("newdata should at least include the following covariates: ","\n",x$Xnames2),collapse=" "))
if (!inherits(newdata, "data.frame")) stop("newdata should be a data.frame object")
#if(plot==TRUE & (missing(mfrow) | missing(x.plot))) stop("Arguments mfrow and x.plot are required if plot=TRUE")

if(x$conv==1|x$conv==2) 
{
  if(x$conv==2 & draws==TRUE)
  {
   cat("No confidence interval will be provided since the program did not converge properly \n")
   draws <- FALSE
  }
  if(!(na.action%in%c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument")
  
  
  ### Traitement des donnees manquantes
  newdata <- newdata[,x$Xnames2]
  newdata <- as.data.frame(newdata)
  colnames(newdata) <- x$Xnames2
  
  linesNA <- apply(newdata,2,function(v) which(is.na(v)))
  linesNA <- unique(unlist(linesNA))
  
  if(length(linesNA) & na.action==2) stop("newdata contains missing values")
  if(length(linesNA) & na.action==1)
  {                                                                              
   newdata <- as.data.frame(newdata[-linesNA,])
   colnames(newdata) <- x$Xnames2[-1]
  }
  
  
  ##pour les facteurs
  olddata <- eval(x$call$data)
  termes <- x$Xnames[-1]
  
   #cas où une variable dans le dataset du modèle est un facteur
   for(v in x$Xnames2[-1])
   {
    if (is.factor(olddata[,v]))
    {
     mod <- levels(olddata[,v])
     if (!(levels(as.factor(newdata[,v])) %in% mod)) stop(paste("invalid level in factor", v))
     newdata[,v] <- factor(newdata[,v], levels=mod)
  
     for(m in mod)
     {
      termes <- gsub(paste(v,m,sep=""),v,termes)
     }
     
    }
   }
  
   #cas où on a factor() dans l'appel de la fonction
   dans_appel <- c(all.names(x$call$fixed),all.names(x$call$random),all.names(x$call$mixture),all.names(x$call$classmb))
   ind_factor <- which(dans_appel=="factor")
   if(length(ind_factor))
   {
    nom.factor <- dans_appel[ind_factor+1]
    for (v in nom.factor)
    {
     mod <- levels(as.factor(olddata[,v]))
     if (!all(levels(as.factor(newdata[,v])) %in% mod)) stop(paste("invalid level in factor", v))
     newdata[,v] <- factor(newdata[,v], levels=mod)
     
     factorv <- paste("factor(",v,")",sep="")
     for(m in mod)
     {
      termes <- gsub(paste(factorv,m,sep=""),v,termes,fixed=TRUE)  
     }
    }
   }
  
  
  ###matrice avec toutes les var et toutes les interactions
  newdata1 <- model.matrix(as.formula(paste("~",paste(termes,collapse="+"))),data=newdata)
   #remettre les termes dans le bon ordre
   Xnames <- x$Xnames[-1]
   z <- grep("factor\\(",Xnames) 
   if (length(z))
   {
    Xnames <- gsub("factor\\(","",Xnames)
    Xnames[z] <- gsub("\\)","",Xnames[z])
   }
   newdata1 <- newdata1[,c("(Intercept)",Xnames)]
  
  
  #newdata1 <- model.matrix(as.formula(paste("~",paste(x$Xnames[-1],collapse="+"))),data=newdata)
  
  #arguments pour predictMult
   maxmes <- dim(newdata1)[1]
   X0 <- as.vector(newdata1)
   ncor <- x$N[7]
   nalea <- x$N[6]
   nv <- dim(newdata1)[2]                          
   ny <- x$N[8]
   best <- x$best
   nef <- x$N[3]
   nvc <- x$N[4]
   if(nvc>0) best[nef+1:nvc] <- x$cholesky[-1]
   npm <- length(best)
   nbzitr <- rep(2,ny)
   nbzitr[which(x$linktype==2)] <- x$nbnodes
   Ymarg <- matrix(0,maxmes*ny,x$ng)
   
   #if(verbose==TRUE) print(head(newdata1))
  
  if(!draws)
  { 
   #if(verbose==TRUE) cat("ny=",ny,"nvc=",nvc,"ncontr=",x$N[2],"nv=",nv,"\n",x$Xnames,"\n",head(newdata1),"\n","idg",x$idg0)
   
   out <- .Fortran("predict_mult",as.double(X0),as.integer(x$idprob0),as.integer(x$idea0),
   as.integer(x$idg0),as.integer(x$idcor0),as.integer(x$idcontr0),as.integer(x$ng),as.integer(ncor),
   as.integer(nalea),as.integer(nv),as.integer(ny),as.integer(maxmes),as.integer(x$idiag),as.integer(x$N[5]),
   as.integer(npm),as.double(best),as.double(x$epsY),as.integer(x$linktype),as.integer(nbzitr),
   as.double(x$linknodes),as.integer(nsim),as.integer(methInteg),Ymarg=as.double(Ymarg),PACKAGE="lcmm")
  
   out$Ymarg[out$Ymarg==9999] <- NA
   #Ypred <- matrix(out$Ymarg,ncol=x$ng,byrow=FALSE)
   Ypred <- data.frame(rep(x$Ynames,each=maxmes),matrix(out$Ymarg,ncol=x$ng,byrow=FALSE))
  
   if (x$ng==1) colnames(Ypred) <- c("Yname","Ypred")
   if (x$ng>1) colnames(Ypred) <- c("Yname",paste("Ypred_class",1:x$ng,sep=""))
    
#   if(plot==TRUE)
#   {
#    if(missing(file)) {dev.new()}
#    if(!missing(file)) {postscript(file,width=6,height=6,paper="special",horizontal=FALSE)}
#    par(mfrow=mfrow, oma=c(0,0,2,0))
#    npred <- dim(newdata)[1]
#    for (i in 1:length(x$Ynames))
#    { plot(newdata[,x.plot], Ypred[1:npred+(i-1)*npred,1],type="l",col=1,xlab=x.plot,ylab=x$Ynames[i],ylim=c(min(Ypred[1:npred+(i-1)*npred,]),max(Ypred[1:npred+(i-1)*npred,])))
#      if(x$ng>1)
#      {
#       for (j in 2:x$ng)
#       {
#        lines(newdata[,x.plot], Ypred[1:npred+(i-1)*npred,j],col=j)
#       }
#       legend("bottomleft",legend=paste("class",1:x$ng),col=1:x$ng,lty=1,bty="n")
#      }
#    }
#    title("Marginal predictions in the outcome scale", outer=TRUE)
#    if(!missing(file)) {dev.off()}  
#   }
  
  }
  else  #draws
  { 
   ndraws <- as.integer(ndraws)
   ydraws <- NULL
   
  
   Mat <- matrix(0,ncol=npm,nrow=npm)
   Mat[upper.tri(Mat,diag=TRUE)]<- x$V
   Chol <- chol(Mat)
   Chol <- t(Chol)
    
   for (j in 1:ndraws)
   {   #cat("boucle sur ndraws j=",j,"\n")
    bdraw <- rnorm(npm)
    bdraw <- x$best + Chol %*% bdraw
   
    out <- .Fortran("predict_mult",as.double(X0),as.integer(x$idprob0),as.integer(x$idea0),
    as.integer(x$idg0),as.integer(x$idcor0),as.integer(x$idcontr0),as.integer(x$ng),as.integer(ncor),
    as.integer(nalea),as.integer(nv),as.integer(ny),as.integer(maxmes),as.integer(x$idiag),as.integer(x$N[5]),
    as.integer(npm),as.double(bdraw),as.double(x$epsY),as.integer(x$linktype),as.integer(nbzitr),
    as.double(x$linknodes),as.integer(nsim),as.integer(methInteg),Ymarg=as.double(Ymarg),PACKAGE="lcmm")
  
    out$Ymarg[out$Ymarg==9999] <- NA
    ydraws <- cbind(ydraws,out$Ymarg)
   }
  
   f <- function(x) {
   quantile(x[!is.na(x)],probs=c(0.025,0.5,0.975))
   }
   ydistr <- apply(ydraws,1,FUN=f)
   Ypred_50 <- matrix(ydistr[2,],ncol=x$ng,byrow=F)
   Ypred_2.5 <- matrix(ydistr[1,],ncol=x$ng,byrow=F)
   Ypred_97.5 <- matrix(ydistr[3,],ncol=x$ng,byrow=F)
  
   #Ypred <- cbind(Ypred_50,Ypred_2.5,Ypred_97.5)
   Ypred <- data.frame(rep(x$Ynames,each=maxmes),Ypred_50,Ypred_2.5,Ypred_97.5)
   
   if (x$ng==1) colnames(Ypred) <- c("Yname","Ypred_50","Ypred_2.5","Ypred_97.5")
   if (x$ng>1) colnames(Ypred) <- c("Yname",c(paste("Ypred_50_class",1:x$ng,sep=""),paste("Ypred_2.5_class",1:x$ng,sep=""),paste("Ypred_97.5_class",1:x$ng,sep="")))
    
#   if(plot==TRUE)
#   {
#    if(missing(file)) {dev.new()}
#    if(!missing(file)) {postscript(file,width=6,height=6,paper="special",horizontal=FALSE)}
#    par(mfrow=mfrow, oma=c(0,0,2,0))
#    npred <- dim(newdata)[1]
#    for (i in 1:length(x$Ynames))
#    { 
#     plot(newdata[,x.plot], Ypred[1:npred+(i-1)*npred,1],type="l",lty=1,col=1,xlab=x.plot,ylab=x$Ynames[i],ylim=c(min(Ypred[1:npred+(i-1)*npred,1:x$ng]),max(Ypred[1:npred+(i-1)*npred,(2*x$ng+1):(3*x$ng)])))
#     lines(newdata[,x.plot], Ypred[1:npred+(i-1)*npred,x$ng+1],lty=2,col=1)
#     lines(newdata[,x.plot], Ypred[1:npred+(i-1)*npred,2*x$ng+1],lty=2,col=1)   
#      if(x$ng>1)
#      {
#       for (j in 2:x$ng)
#       {
#        lines(newdata[,x.plot], Ypred[1:npred+(i-1)*npred,j],col=j,lty=1)
#        lines(newdata[,x.plot], Ypred[1:npred+(i-1)*npred,x$ng+j],col=j,lty=2)
#        lines(newdata[,x.plot], Ypred[1:npred+(i-1)*npred,2*x$ng+j],col=j,lty=2)
#       }
#       legend("bottomleft",legend=paste("class",1:x$ng),col=1:x$ng,lty=1,bty="n")
#      }
#    }
#    title("Marginal predictions in the outcome scale", outer=TRUE)  
#    if(!missing(file)) {dev.off()}  
#   }
  }

}

else
{
 cat(" The program stopped abnormally. No prediction can be computed.\n")
 Ypred <- NA
}

 Ypred
}
