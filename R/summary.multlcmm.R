summary.multlcmm <- function(object,...)
{
 x <- object
 if (!inherits(x, "multlcmm")) stop("use only with \"multlcmm\" objects")

 cat("General latent class mixed model", "\n")
 cat("     fitted by maximum likelihood method", "\n")

 cl <- x$call
 cl$B <- NULL
 cat(" \n")
 dput(cl)
 cat(" \n")

 cat("Statistical Model:", "\n")
 cat(paste("     Dataset:", x$call$data),"\n")
 cat(paste("     Number of subjects:", x$ns),"\n")

 cat(paste("     Number of observations:", x$N[9]),"\n")
 cat(paste("     Number of latent classes:", x$ng), "\n")
 cat(paste("     Number of parameters:", length(x$best))," \n")

 ntrtot <- rep(NA,x$N[8])
 numSPL <- 0
 cat("     Link functions: ")
 for (yk in 1:x$N[8])
 {
  if (x$linktype[yk]==0)
  {
   ntrtot[yk] <- 2
   if (yk>1) cat("                     ")
   cat("Linear for",x$Ynames[yk]," \n")
  }
  if (x$linktype[yk]==1)
  {
  ntrtot[yk] <- 4
  if (yk>1) cat("                     ")
  cat("Standardised Beta CdF for",x$Ynames[yk]," \n")
  }
  if (x$linktype[yk]==2) 
  {
   numSPL <- numSPL+1
   ntrtot[yk] <- x$nbnodes[numSPL]+2
   if (yk>1) cat("                     ")
   cat("Quadratic I-splines with nodes", x$linknodes[1:x$nbnodes[numSPL],yk],"for",x$Ynames[yk], "\n")
  }
 }

 cat(" \n")
 cat("Iteration process:", "\n")

 if(x$conv==1) cat("     Convergence criteria satisfied")
 if(x$conv==2) cat("     Maximum number of iteration reached without convergence")
 if(x$conv==4|x$conv==12)
 {
  cat("     The program stopped abnormally. No results can be displayed.\n")
 }
 else
 {
  cat(" \n")
  cat("     Number of iterations: ", x$niter, "\n")
  cat("     Convergence criteria: parameters=", signif(x$gconv[1],2), "\n")
  cat("                         : likelihood=", signif(x$gconv[2],2), "\n")
  cat("                         : second derivatives=", signif(x$gconv[3],2), "\n")
  cat(" \n")
  cat("Goodness-of-fit statistics:", "\n")
  cat(paste("     maximum log-likelihood:", round(x$loglik,2))," \n")
  cat(paste("     AIC:", round(-2*x$loglik+2*length(x$best),2))," \n")
  cat(paste("     BIC:", round(-2*x$loglik+length(x$best)*log(x$ns),2))," \n")
  cat(" \n")


  cat("Maximum Likelihood Estimates:", "\n")
  cat(" \n")
  
  nprob <- x$N[1]
  ncontr   <- x$N[2]
  nef   <- x$N[3]
  nvc    <- x$N[4]
  nw <- x$N[5]
  nalea <- x$N[6]
  ncor <- x$N[7]
  ny <- x$N[8]
  NPM   <- length(x$best)
  
  
  se <- rep(NA,NPM)
  if (x$conv==1)
  {
   #recuperation des indices de V
   id <- 1:NPM
   indice <- id*(id+1)/2
   se <-sqrt(x$V[indice])
   if (nvc>0) se[(nef+1):(nef+nvc)]<-NA
   wald <- x$best/se
   pwald <- 1-pchisq(wald**2,1)
   coef <- x$best
  }
  else
  {
   se <- NA
   wald <- NA
   pwald <- NA
   coef <- x$best
  }
  
  if(nprob>0)
  {
   cat("Fixed effects in the class-membership model:\n" )
   cat("(the class of reference is the last class) \n")
   if (x$ng<3)
   {
    tmp <- cbind(round(coef[1:nprob],5),round(se[1:nprob],5),round(wald[1:nprob],3),round(pwald[1:nprob],5))
    dimnames(tmp) <- list(names(coef)[1:nprob], c("coef", "Se", "Wald", "p-value"))
   }
   else
   {
    tmp <- cbind(round(coef[1:(x$ng-1)],5),round(se[1:(x$ng-1)],5),round(wald[1:(x$ng-1)],3),round(pwald[1:(x$ng-1)],5))
    dimnames(tmp) <- list(names(coef)[1:(x$ng-1)], c("coef", "Se", "Wald", "p-value"))
   
    indice2 <- 1:NPM*(1:NPM+1)/2
    nom.prob <- x$Xnames[as.logical(x$idprob0)]
    for (i in 2:sum(x$idprob0))
    {
     #matrice de variance pour test global
     indtmp <- indice2[((i-1)*(x$ng-1)+1):(i*(x$ng-1))]
     indtmp <- cbind(indtmp-0:(length(indtmp)-1),indtmp)
     indV <- NULL
     for (j in 1:dim(indtmp)[1])
     {
      indV <- c(indV,seq(indtmp[j,1],indtmp[j,2]))
     }
     Vprob <- matrix(0,x$ng-1,x$ng-1)
     Vprob[upper.tri(Vprob,diag=TRUE)] <- x$V[indV]
     Vprob <- t(Vprob)
     Vprob[upper.tri(Vprob)] <- Vprob[lower.tri(Vprob)]
    
     vect.prob <- coef[((i-1)*(x$ng-1)+1):(i*(x$ng-1))]
     wald.prob <- t(vect.prob) %*% solve(Vprob,vect.prob)
     p.wald.prob <- 1-pchisq(wald.prob,x$ng-1)
  
     tmp2 <- cbind(vect.prob,se[((i-1)*(x$ng-1)+1):(i*(x$ng-1))],wald[((i-1)*(x$ng-1)+1):(i*(x$ng-1))],pwald[((i-1)*(x$ng-1)+1):(i*(x$ng-1))])
     tmp2 <- rbind(rep(NA,4),tmp2)
     
     if(x$conv==1)
     {
      z <- paste(" class",1:(x$ng-1))
      rownames(tmp2) <- c(paste(nom.prob[i]," (p=",round(p.wald.prob,5),"):",sep=""),z)
     }
     if(x$conv==2)
     {
      z <- paste(" class",1:(x$ng-1))
      rownames(tmp2) <- c(paste(nom.prob[i],":",sep=""),z)
     }      
    }
    tmp <- rbind(tmp,tmp2)
   }
   cat("\n")
   prmatrix(tmp,na.print="")
   cat("\n")
  }

  cat("Fixed effects in the longitudinal model:\n" )
  
  tmp <- matrix(c(0,NA,NA,NA),nrow=1,ncol=4)
  if (nef>0)
  {
   tmp2 <- cbind(round(coef[(nprob+1):(nef-ncontr)],5),round(se[(nprob+1):(nef-ncontr)],5),round(wald[(nprob+1):(nef-ncontr)],3),round(pwald[(nprob+1):(nef-ncontr)],5))
   tmp <- rbind(tmp,tmp2)
  } 
  interc <- "intercept"
  if (x$ng>1)
  {
   interc <- paste(interc,"class1")
  }
  interc <- paste(interc,"(not estimated)")
  if(nef>0) dimnames(tmp) <- list(c(interc,names(coef)[(nprob+1):(nef-ncontr)]), c("coef", "Se", "Wald", "p-value"))
  else dimnames(tmp) <- list(interc, c("coef", "Se", "Wald", "p-value"))
  cat("\n")
  
  if(ncontr>0)
  {
   indice2 <- 1:NPM*(1:NPM+1)/2
   nom.contr <- x$Xnames[as.logical(x$idcontr0)]
   for (i in 1:sum(x$idcontr0))
   {
   #matrice de variance pour test et se du dernier coef
    indtmp <- indice2[(nef-ncontr+(i-1)*(ny-1)+1):(nef-ncontr+i*(ny-1))]
    indtmp <- cbind(indtmp-0:(length(indtmp)-1),indtmp)
    indV <- NULL
    for (j in 1:dim(indtmp)[1])
    {
     indV <- c(indV,seq(indtmp[j,1],indtmp[j,2]))
    }
    Vcontr <- matrix(0,ny-1,ny-1)
    Vcontr[upper.tri(Vcontr,diag=TRUE)] <- x$V[indV]
    Vcontr <- t(Vcontr)
    Vcontr[upper.tri(Vcontr)] <- Vcontr[lower.tri(Vcontr)]
    
    vect.gamma <- coef[(nef-ncontr+(i-1)*(ny-1)+1):(nef-ncontr+i*(ny-1))]
    wald.contr <- t(vect.gamma) %*% solve(Vcontr,vect.gamma)
    p.wald.contr <- 1-pchisq(wald.contr,ny-1)
  
    tmp2 <- cbind(vect.gamma,se[(nef-ncontr+(i-1)*(ny-1)+1):(nef-ncontr+i*(ny-1))],wald[(nef-ncontr+(i-1)*(ny-1)+1):(nef-ncontr+i*(ny-1))],pwald[(nef-ncontr+(i-1)*(ny-1)+1):(nef-ncontr+i*(ny-1))])
    tmp2 <- rbind(rep(NA,4),tmp2)
  
    if(x$conv==1)
    {
     pp <- -sum(na.omit(tmp2[,1]))/sqrt(sum(Vcontr))
     tmp2 <- rbind(tmp2,c(-sum(na.omit(tmp2[,1])),sqrt(sum(Vcontr)),pp,1-pchisq(pp*pp,1)))
     if(round(p.wald.contr,5)!=0)rownames(tmp2) <- c(paste("Contrasts on ",nom.contr[i]," (p=",round(p.wald.contr,5),")",sep=""),x$Ynames)
     if(round(p.wald.contr,5)==0)rownames(tmp2) <- c(paste("Contrasts on ",nom.contr[i]," (p<0.00001)",sep=""),x$Ynames)
    }
    if(x$conv==2)
    {
     tmp2 <- rbind(tmp2,c(-sum(na.omit(tmp2[,1])),NA,NA,NA))
     rownames(tmp2) <- c(paste("Contrasts on ",nom.contr[i],sep=""),x$Ynames)
    }
    rownames(tmp2)[nrow(tmp2)] <- paste(rownames(tmp2)[nrow(tmp2)],"*",sep="")
    tmp <- rbind(tmp,tmp2)
   }
  }
  
  #res <- tmp
  prmatrix(round(tmp,5),na.print="")
  cat("\n")
  if(ncontr>0)
  {
   cat(" *coefficient not estimated but obtained from the others as minus the sum of them")
   cat("\n")
  }
  
  cat("\n")
  cat("Variance-covariance matrix of the random-effects:\n" )
  cat("(the variance of the first random effect is not estimated)\n")
  if(x$idiag==1)
  {
   if (nvc>0) 
   {
    Mat.cov <- diag(c(1,coef[(nef+1):(nef+nvc)]))
   }
   else
   {
    Mat.cov <- matrix(1,ncol=1)
   }
   Mat.cov[lower.tri(Mat.cov)] <- 0
   Mat.cov[upper.tri(Mat.cov)] <- NA
  }
  if(x$idiag==0)
  {
   Mat.cov<-matrix(0,ncol=sum(x$idea0),nrow=sum(x$idea0))
   if(nvc>0) 
   {
    Mat.cov[upper.tri(Mat.cov,diag=TRUE)]<-c(1,coef[(nef+1):(nef+nvc)])
    Mat.cov <-t(Mat.cov)
    Mat.cov[upper.tri(Mat.cov)] <- NA
   }
   else Mat.cov[1,1] <- 1
  }
  colnames(Mat.cov) <-x$Xnames[x$idea0==1]
  rownames(Mat.cov) <-x$Xnames[x$idea0==1]
  print(round(Mat.cov,5),na.print="")
  cat("\n")
  
  std <- NULL
  nom <- NULL
  if(nw>=1) 
  {
   nom <- paste("Proportional coefficient class",c(1:(x$ng-1)),sep="")
   std <-cbind(coef[(nef+nvc+1):(nef+nvc+nw)],se[(nef+nvc+1):(nef+nvc+nw)])
  }
  if(ncor==2)
  {
   nom <- c(nom,"AR correlation parameter:","AR standard error:")
   std <-rbind(std,c(coef[(nef+nvc+nw+1)],se[(nef+nvc+nw+1)]),c(abs(coef[(nef+nvc+nw+2)]),se[(nef+nvc+nw+2)]))
  }
  if(ncor==1) 
  {
   nom <- c(nom,"BM standard error:")
   std <-rbind(std,c(abs(coef[(nef+nvc+nw+1)]),se[(nef+nvc+nw+1)]))
  }
  if (!is.null(std)) 
  {
   rownames(std) <- nom
   colnames(std) <-c("coef","Se")
   print(round(std,5), na.print="")
   cat("\n")
  }

  
  std.err <- NULL
  nom <- NULL
  std.err <- rbind(std.err,abs(x$best[nef+nvc+nw+ncor+1:ny]))
  nom <- c(nom, "Residual standard error:")
  if(nalea>0)
  {
   std.err <- rbind(std.err,abs(x$best[nef+nvc+nw+ncor+ny+1:nalea]))
   nom <- c(nom, "Standard error of the random effect:")
  }
  colnames(std.err) <- x$Ynames
  rownames(std.err) <- nom
  print(round(std.err,5))
  
  cat("\n")
  
  cat("Parameters of the link functions:\n" )
  
  tmp <- cbind(round(coef[(nef+nvc+nw+ncor+ny+nalea+1):NPM],5),round(se[(nef+nvc+nw+ncor+ny+nalea+1):NPM],5),round(wald[(nef+nvc+nw+ncor+ny+nalea+1):NPM],3),round(pwald[(nef+nvc+nw+ncor+ny+nalea+1):NPM],5))
  colnames(tmp) <- c("coef", "Se", "Wald", "p-value")
  tmp.rownames <- NULL
  for (yk in 1:ny)
  {
   tmp.rownames <- c(tmp.rownames, paste(rep(x$Ynames[yk],ntrtot[yk]),names(coef[(nef+nvc+nw+ncor+ny+nalea+sum(ntrtot[1:yk])-ntrtot[yk]+1):(nef+nvc+nw+ncor+ny+nalea+sum(ntrtot[1:yk]))]),sep="-"))
  }
  rownames(tmp) <- tmp.rownames
  cat("\n")
  prmatrix(tmp,na.print="")
  cat("\n")
 }
 #invisible(res)
}
