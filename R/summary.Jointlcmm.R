summary.Jointlcmm <- function(object,...)
{
 x <- object
 if (!inherits(x, "Jointlcmm")) stop("use only with \"Jointlcmm\" objects")

 cat("Joint latent class model for quantitative outcome and competing risks", "\n")
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

  nbevt <- length(x$hazard[[1]])
  nprisq <- rep(NA,nbevt)
  nrisq <- rep(NA,nbevt)

  typrisq <- x$hazard[[1]]
  hazardtype <- x$hazard[[2]]
  nz <- x$hazard[[4]] 

 for(ke in 1:nbevt)
     {
        if(typrisq[ke]==1) nprisq[ke] <- nz[ke]-1
        if(typrisq[ke]==2) nprisq[ke] <- 2
        if(typrisq[ke]==3) nprisq[ke] <- nz[ke]+2

        if(hazardtype[ke]=="Common") nrisq[ke] <- nprisq[ke]
        if(hazardtype[ke]=="PH") nrisq[ke] <- nprisq[ke]+x$ng-1
        if(hazardtype[ke]=="Specific") nrisq[ke] <- nprisq[ke]*x$ng

        cat(paste("     Event ",ke,": \n",sep=""))
        cat(paste("        Number of events: ", x$N[9+ke],"\n",sep=""))
        if(x$ng>1)
            {
                if (hazardtype[ke]=="Specific") cat("        Class-specific hazards and \n")
                if (hazardtype[ke]=="PH") cat("        Proportional hazards over latent classes and \n")
                if (hazardtype[ke]=="Common") cat("        Common hazards over classes and \n")
            }
        
        if (typrisq[ke]==2)
            {
                cat("        Weibull baseline risk function \n")
            }
        if (typrisq[ke]==1)
            {
                cat("        Piecewise constant baseline risk function with nodes \n")
                cat("       ",x$hazard[[3]][1:nz[ke],ke]," \n")
            }
        if (typrisq[ke]==3)
            {
                cat("        M-splines constant baseline risk function with nodes \n")
                cat("       ",x$hazard[[3]][1:nz[ke],ke]," \n")
            }
        
        
    }

 ntrtot <- x$N[8]
 numSPL <- 0
 if(x$linktype!=-1)
     {
         cat(paste("     Link function for ",x$Names$Yname,": ",sep=""))
         if (x$linktype==0)
             {
                 cat("Linear \n")
             }
         if (x$linktype==1)
             {
                 cat("Standardised Beta CdF \n")
             }
         if (x$linktype==2) 
             {
                 cat("Quadratic I-splines with nodes ", x$linknodes ,"\n")
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
  if(!is.na(x$scoretest[1])&(length(x$hazard[[1]])==1)){
    cat(paste("     Score test statistic for CI assumption: ", round(x$scoretest[1],3)," (p-value=",round((1-pchisq(x$scoretest[1],sum(x$idea))),4),")" ,sep=""))
  }
  if(!is.na(x$scoretest[1])&(length(x$hazard[[1]])>1)){
    cat(paste("     Score test statistic for global CI assumption: ", round(x$scoretest[1],3)," (p-value=",round((1-pchisq(x$scoretest[1],sum(x$idea))),4),")" ,sep=""),"\n")
  }
  if(!is.na(x$scoretest[1])&(length(x$hazard[[1]])>1)){
    cat("     Score test statistic for event-specific CI assumption: \n")
       for (ke in 1:length(x$hazard[[1]])){ 
        if(!is.na(x$scoretest[1+ke])){
        cat(paste("           event ",ke,":", round(x$scoretest[1+ke],3)," (p-value=",round((1-pchisq(x$scoretest[1+ke],sum(x$idea))),4),")" ,sep=""),"\n")
  }
  else{
    cat(paste("           event ",ke,": problem in the computation", "\n"))
  }
            }
  }
  
  cat(" \n")
  cat(" \n")

  cat("Maximum Likelihood Estimates:", "\n")
  cat(" \n")
  
  nprob <- x$N[1]
  nrisqtot <- x$N[2]
  nvarxevt <- x$N[3]
  nef <- x$N[4]
  nvc <- x$N[5]
  nw <- x$N[6]
  ncor <- x$N[7]
  ntrtot <- x$N[8]
  NPM <- length(x$best)

  #nvdepsurv <- length(x$Name$TimeDepVar.name)

      
  se <- rep(NA,NPM)
  if (x$conv==1)
  {
   #recuperation des indices de V
   id <- 1:NPM
   indice <- id*(id+1)/2
   se <-sqrt(x$V[indice])
   #if (nvc>0) se[(nef+1):(nef+nvc)]<-NA ??
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
  ## if (x$ng<3)
  ## {
    tmp <- cbind(round(coef[1:nprob],5),round(se[1:nprob],5),
                 round(wald[1:nprob],3),round(pwald[1:nprob],5))
    dimnames(tmp) <- list(names(coef)[1:nprob], c("coef", "Se", "Wald", "p-value"))
   ## }
   ## else
   ## {
   ##  tmp <- cbind(round(coef[1:(x$ng-1)],5),round(se[1:(x$ng-1)],5),
   ##               round(wald[1:(x$ng-1)],3),round(pwald[1:(x$ng-1)],5))
   ##  dimnames(tmp) <- list(names(coef)[1:(x$ng-1)], c("coef", "Se", "Wald", "p-value"))
   
   ##  indice2 <- 1:NPM*(1:NPM+1)/2
   ##  nom.prob <- x$Names$Xnames[as.logical(x$idprob)]
   ##  for (i in 2:sum(x$idprob))
   ##  {
   ##   #matrice de variance pour test global
   ##   indtmp <- indice2[((i-1)*(x$ng-1)+1):(i*(x$ng-1))]
   ##   indtmp <- cbind(indtmp-0:(length(indtmp)-1),indtmp)
   ##   indV <- NULL
   ##   for (j in 1:dim(indtmp)[1])
   ##   {
   ##    indV <- c(indV,seq(indtmp[j,1],indtmp[j,2]))
   ##   }
   ##   Vprob <- matrix(0,x$ng-1,x$ng-1)
   ##   Vprob[upper.tri(Vprob,diag=TRUE)] <- x$V[indV]
   ##   Vprob <- t(Vprob)
   ##   Vprob[upper.tri(Vprob)] <- Vprob[lower.tri(Vprob)]
    
   ##   vect.prob <- coef[((i-1)*(x$ng-1)+1):(i*(x$ng-1))]
   ##   wald.prob <- t(vect.prob) %*% solve(Vprob,vect.prob)
   ##   p.wald.prob <- 1-pchisq(wald.prob,x$ng-1)
  
   ##   tmp2 <- cbind(vect.prob,se[((i-1)*(x$ng-1)+1):(i*(x$ng-1))],
   ##                 wald[((i-1)*(x$ng-1)+1):(i*(x$ng-1))],
   ##                 pwald[((i-1)*(x$ng-1)+1):(i*(x$ng-1))])
   ##   tmp2 <- rbind(rep(NA,4),tmp2)
     
   ##   if(x$conv==1)
   ##   {
   ##    z <- paste(" class",1:(x$ng-1))
   ##    rownames(tmp2) <- c(paste(nom.prob[i]," (p=",round(p.wald.prob,5),"):",sep=""),z)
   ##   }
   ##   if(x$conv==2)
   ##   {
   ##    z <- paste(" class",1:(x$ng-1))
   ##    rownames(tmp2) <- c(paste(nom.prob[i],":",sep=""),z)
   ##   }      
   ##  }
   ##  tmp <- rbind(tmp,tmp2)
   ## }
   cat("\n")
   prmatrix(tmp,na.print="")
   cat("\n")
}


cat("Parameters in the proportional hazard model:\n" )

tmp <- cbind(round(coef[nprob+1:(nrisqtot+nvarxevt)],5),
             round(se[nprob+1:(nrisqtot+nvarxevt)],5),
             round(wald[nprob+1:(nrisqtot+nvarxevt)],3),
             round(pwald[nprob+1:(nrisqtot+nvarxevt)],5))
dimnames(tmp) <- list(names(coef)[nprob+1:(nrisqtot+nvarxevt)],
                      c("coef", "Se", "Wald", "p-value"))
cat("\n")
prmatrix(tmp,na.print="")
cat("\n")

  

  cat("Fixed effects in the longitudinal model:\n" )
  
  if(x$linktype!=-1) tmp <- matrix(c(0,NA,NA,NA),nrow=1,ncol=4)
  if(x$linktype==-1) tmp <- NULL
  if (nef>0)
  {
   tmp2 <- cbind(round(coef[nprob+nrisqtot+nvarxevt+1:nef],5),
                 round(se[nprob+nrisqtot+nvarxevt+1:nef],5),
                 round(wald[nprob+nrisqtot+nvarxevt+1:nef],3),
                 round(pwald[nprob+nrisqtot+nvarxevt+1:nef],5))
   tmp <- rbind(tmp,tmp2)
  }
  
  interc <- "intercept"
  if (x$ng>1)
  {
   interc <- paste(interc,"class1")
  }
  if(x$linktype!=-1) interc <- paste(interc,"(not estimated)")
  if(x$linktype==-1) interc <- NULL
  
  if(nef>0) dimnames(tmp) <- list(c(interc,names(coef)[nprob+nrisqtot+nvarxevt+1:nef]), c("coef", "Se", "Wald", "p-value"))
  else dimnames(tmp) <- list(interc, c("coef", "Se", "Wald", "p-value"))
  cat("\n")

  tTable <- tmp
  
  prmatrix(round(tmp,5),na.print="")
  cat("\n")

  if(nvc>0)
  {    
   cat("\n")
   cat("Variance-covariance matrix of the random-effects:\n" )
   if(x$idiag==1)
   {
    Mat.cov <- diag(coef[nprob+nrisqtot+nvarxevt+nef+1:nvc])
    Mat.cov[lower.tri(Mat.cov)] <- 0
    Mat.cov[upper.tri(Mat.cov)] <- NA
    if(nvc==1) Mat.cov <- matrix(coef[nprob+nrisqtot+nvarxevt+nef+1:nvc],1,1)
   }
   if(x$idiag==0)
   {
    Mat.cov<-matrix(0,ncol=sum(x$idea),nrow=sum(x$idea))
    Mat.cov[upper.tri(Mat.cov,diag=TRUE)] <- coef[nprob+nrisqtot+nvarxevt+nef+1:nvc]
    Mat.cov <-t(Mat.cov)
    Mat.cov[upper.tri(Mat.cov)] <- NA
   }
 
   colnames(Mat.cov) <-x$Names$Xnames[x$idea==1]
   rownames(Mat.cov) <-x$Names$Xnames[x$idea==1]
   prmatrix(round(Mat.cov,5),na.print="")
   cat("\n")
  }


  
  std <- NULL
  nom <- NULL
  if(nw>=1) 
  {
   nom <- paste("Proportional coefficient class",c(1:(x$ng-1)),sep="")
   std <-cbind(coef[nprob+nrisqtot+nvarxevt+nef+nvc+1:nw],
               se[nprob+nrisqtot+nvarxevt+nef+nvc+1:nw])
  }
  if(ncor==2)
  {
   nom <- c(nom,"AR correlation parameter:","AR standard error:")
   std <-rbind(std,c(coef[nprob+nrisqtot+nvarxevt+nef+nvc+nw+1],
                     se[nprob+nrisqtot+nvarxevt+nef+nvc+nw+2]),
                   c(abs(coef[nprob+nrisqtot+nvarxevt+nef+nvc+nw+2]),
                     se[nprob+nrisqtot+nvarxevt+nef+nvc+nw+2]))
  }
  if(ncor==1) 
  {
   nom <- c(nom,"BM standard error:")
   std <-rbind(std,c(abs(coef[nprob+nrisqtot+nvarxevt+nef+nvc+nw+1]),
                     se[nprob+nrisqtot+nvarxevt+nef+nvc+nw+1]))
  }
  if (!is.null(std)) 
  {
   rownames(std) <- nom
   colnames(std) <- c("coef","Se")
   prmatrix(round(std,5), na.print="")
   cat("\n")
  }

  if(x$linktype==-1)
      {
          tmp <- cbind(round(coef[NPM],5),round(se[NPM],5))
          rownames(tmp) <- "Residual standard error"
          colnames(tmp) <- c("coef","se")
          prmatrix(tmp, na.print="")
          cat("\n")
      }
  else
      {
          cat("Residual standard error (not estimated) = 1\n")
          cat("\n")
          
  
          cat("Parameters of the link function:\n" )
  
          tmp <- cbind(round(coef[(nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor+1):NPM],5),
                       round(se[(nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor+1):NPM],5),
                       round(wald[(nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor+1):NPM],3),
                       round(pwald[(nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor+1):NPM],5))
          colnames(tmp) <- c("coef", "Se", "Wald", "p-value")
          rownames(tmp) <- names(x$best[(nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor+1):NPM])
          cat("\n")
          prmatrix(tmp,na.print="")
          cat("\n")
      }
  
  return(invisible(tTable))
 }
}
