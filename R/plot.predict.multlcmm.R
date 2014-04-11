plot.predict.multlcmm <- function(x,newdata,var.time,legend.loc="topright",na.action=1,legend,add=FALSE,confint=FALSE,...)
{
 if(missing(var.time)) stop("The argument var.time should be specified")
 if(missing(newdata)) stop("The argument newdata should be specified")
 if(missing(x)) stop("The argument x should be specified")
 if (!inherits(x, "multlcmm")) stop("use only with \"multlcmm\" objects")
 if (!all(x$Xnames2 %in% colnames(newdata))) {
 stop(paste(c("newdata should at least include the following covariates: ","\n",x$Xnames2),collapse=" "))}
 if (!inherits(newdata, "data.frame")) stop("newdata should be a data.frame object")
 if (!inherits(var.time, "character")) stop("the class of var.time should be character")
 if (!(var.time %in% colnames(newdata))) stop("The variable names var.time should be included in newdata")
 if(is.na(as.logical(add))) stop("add should be TRUE or FALSE")
 
 if(x$conv==1|x$conv==2)
 {
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
  
  ###calcul des predictions
     
  X1 <- NULL
  X2 <- NULL
  b1 <- NULL
  b2 <- NULL
  
  placeV <- list() #places pour les variances
  placeV$commun <- NA
  for(i in 1:x$ng)
  {
   placeV[paste("class",i,sep="")] <- NA
  }  
  
  kk<-0
  for(k in 1:length(x$idg0))
  {
   if(x$idg0[k]==1)
   {
    X1 <- cbind(X1,newdata1[,k])
    if (k==1) b1 <- c(b1,0)
    if (k>1) 
    {
     place <- x$N[1]+kk
     b1 <- c(b1,x$best[place+1])
     placeV$commun <- c(placeV$commun,place+1)      
     kk <- kk+1
    }
   }
  
   if(x$idg0[k]==2)
   {
    X2 <- cbind(X2,newdata1[,k])
    if (k==1)
    {
     place1 <- x$N[1]+kk+1
     place2 <- x$N[1]+kk+x$ng-1
     b2 <- rbind(b2,c(0,x$best[place1:place2]))
     for(i in 2:x$ng)
     {
      placeV[[paste("class",i,sep="")]] <- c(placeV[[paste("class",i,sep="")]],x$N[1]+kk+i-1)
     }      
     kk <- kk+x$ng-1
    }
    if (k>1)
    {
     place1 <- x$N[1]+kk+1
     place2 <- x$N[1]+kk+x$ng
     b2 <- rbind(b2,x$best[place1:place2])
     for(i in 1:x$ng)
     {
      placeV[[paste("class",i,sep="")]] <- c(placeV[[paste("class",i,sep="")]],x$N[1]+kk+i)
     }     
     kk <- kk+x$ng
    }
   }
  }
  
  Y<-matrix(0,length(newdata1[,1]),x$ng)
  colnames(Y) <- paste("class",1:x$ng,sep="")
  for(g in 1:x$ng)
  {
   if(length(b1) != 0)
   {
    Y[,g]<- X1 %*% b1
   }
   if(length(b2) != 0)
   {
    Y[,g]<- Y[,g] + X2 %*% b2[,g]
   }
  }

  X <- newdata1[,var.time]
  
  ny <- length(x$Ynames)

  #extraction de Var(beta) 
  Vbeta <- matrix(0,x$N[3]-x$N[2]-x$N[1],x$N[3]-x$N[2]-x$N[1])
  npm <- length(x$best)
  indice <- 1:npm * (1:npm+1) /2
  indtmp <- indice[(x$N[1]+1):(x$N[3]-x$N[2])]
  indtmp <- cbind(indtmp-0:(length(indtmp)-1),indtmp)
  
  indV <- NULL
  for(i in 1:nrow(indtmp))
  {
   indV <- c(indV,seq(indtmp[i,1],indtmp[i,2]))
  }
  
  Vbeta[upper.tri(Vbeta, diag=TRUE)] <- x$V[indV]
  Vbeta <- t(Vbeta)
  Vbeta[upper.tri(Vbeta, diag=TRUE)] <- x$V[indV] #c'est Var(beta)
  
  
  #IC pour les predictions 
  lower <- matrix(0,nrow(Y),ncol(Y))  
  upper <- matrix(0,nrow(Y),ncol(Y))
  colnames(lower) <- paste("lower.class",1:x$ng,sep="")
  colnames(upper) <- paste("upper.class",1:x$ng,sep="") 
   
  if(x$ng==1)
  { 
   #varpred <- diag(X1[,-1,drop=FALSE] %*% Vbeta %*% t(X1[,-1,drop=FALSE]))
   varpred <- apply(X1[,-1,drop=FALSE],1,function(x) matrix(x,nrow=1) %*% Vbeta %*% matrix(x,ncol=1)) 
     # browser()
   lower[,1] <- Y[,1] -1.96 * sqrt(varpred)
   upper[,1] <- Y[,1] +1.96 * sqrt(varpred)
  }
  else
  {
   for(g in 1:x$ng)
   {
    ind <- na.omit(c(placeV[["commun"]],placeV[[paste("class",g,sep="")]]))

    if(g==1)
    {
     if(x$idg0[1]==1)
     {
      X12 <- X12 <- cbind(X1[,-1,drop=FALSE],X2)
     }
     
     if(x$idg0[1]==2)
     {
      X12 <- X12 <- cbind(X1,X2[,-1,drop=FALSE])
     }
    }
    else
    {
     X12 <- cbind(X1,X2)    
    }
    
    X12 <- X12[,order(ind),drop=FALSE]

    varclass <- Vbeta[sort(ind)-x$N[1],sort(ind)-x$N[1]]
    #varpred <- diag(X12 %*% varclass %*% t(X12))
    varpred <- apply(X12,1,function(x) matrix(x,nrow=1) %*% varclass %*% matrix(x,ncol=1)) 
    
    lower[,g] <- Y[,g] -1.96 * sqrt(varpred)
    upper[,g] <- Y[,g] +1.96 * sqrt(varpred)    
   }    
  }  


  dots <- list(...)
  dots <- dots[setdiff(names(dots),c("x","y","log"))]

  if(length(list(...)$main))
  {
   title1 <- as.character(eval(match.call()$main))
   dots <- dots[setdiff(names(dots),"main")]
  }
  else title1 <- "Class-specific mean latent process predicted trajectory"

  if(length(list(...)$col))
  {
   color <- as.vector(eval(match.call()$col))
   dots <- dots[-which(names(dots)=="col")]
  }
  else  color <- 1:x$ng

  if(length(list(...)$type))
  {
   type1 <- eval(match.call()$type)
   dots <- dots[-which(names(dots)=="type")]
  }
  else  type1 <- "l"

  if(length(list(...)$lty))
  {
   lty1 <- eval(match.call()$lty)
   dots <- dots[-which(names(dots)=="lty")]
  }
  else  lty1 <- 1

  if(length(list(...)$ylab))
  {
   ylab1 <- as.character(eval(match.call()$ylab))
   dots <- dots[setdiff(names(dots),"ylab")]
  }
  else ylab1 <- "outcome"

  if(length(list(...)$xlab))
  {
   xlab1 <- as.character(eval(match.call()$xlab))
   dots <- dots[setdiff(names(dots),"xlab")]
  }
  else xlab1 <- var.time
  
  if(missing(legend)) legend <- paste("class",1:x$ng,sep="")

  if(length(list(...)$box.lty))
  {
   box.lty1 <- as.character(eval(match.call()$box.lty))
   dots <- dots[setdiff(names(dots),"box.lty")]
  }
  else box.lty1 <- 0

  if(length(list(...)$inset))
  {
   inset1 <- eval(match.call()$inset)
   dots <- dots[setdiff(names(dots),"inset")]
  }
  else inset1 <- c(0.02,0.02)

  if(isTRUE(confint))
  {
   if(length(list(...)$ylim))
   {
    ylim1 <- eval(match.call()$ylim)
    dots <- dots[setdiff(names(dots),"ylim")]
   }
   else
   {
    ylim1 <- range(c(Y,lower,upper),na.rm=TRUE)
   }
  }
  else
  {
   if(length(list(...)$ylim))
   {
    ylim1 <- eval(match.call()$ylim)
    dots <- dots[setdiff(names(dots),"ylim")]
   }
   else
   {
    ylim1 <- range(Y,na.rm=TRUE)
   }  
  }
  
  names.plot <- c("adj","ann","asp","axes","bg","bty","cex","cex.axis","cex.lab","cex.main","cex.sub","col","col.axis",
  "col.lab","col.main","col.sub","crt","err","family","fig","fin","font","font.axis","font.lab","font.main","font.sub",
  "frame.plot","lab","las","lend","lheight","ljoin","lmitre","lty","lwd","mai","main","mar","mex","mgp","mkh","oma",
  "omd","omi","pch","pin","plt","ps","pty","smo","srt","sub","tck","tcl","type","usr","xaxp","xaxs","xaxt","xlab",
  "xlim","xpd","yaxp","yaxs","yaxt","ylab","ylbias","ylim") 
  dots.plot <- dots[intersect(names(dots),names.plot)]
  
  if(!isTRUE(add))
  {
   do.call("matplot",c(dots.plot,list(x=X,y=Y,xlab=xlab1,ylab=ylab1,main=title1,type=type1,col=color,lty=lty1,ylim=ylim1)))
  }
  else
  {
   do.call("matlines",c(dots.plot,list(x=X,y=Y,type=type1,col=color,lty=lty1))) 
  }
  
  if(isTRUE(confint))
  {
   if(length(list(...)$lwd))
   {
    lwd1 <- eval(match.call()$lwd)
   }
   else
   {
    lwd1 <- 1
   }
        
   matlines(x=X,y=lower,col=color,lwd=lwd1,lty=2)
   matlines(x=X,y=upper,col=color,lwd=lwd1,lty=2)
  }    
   
  names.legend <- c("fill","border","lty","lwd","pch","angle","density","bg","box.lwd",   
  "box.lty","box.col","pt.bg","cex","pt.cex","pt.lwd","xjust","yjust","x.intersp","y.intersp","adj","text.width",
  "text.col","text.font","merge","trace","plot","ncol","horiz","title","xpd","title.col","title.adj","seg.len")  
    
  dots.leg <- dots[intersect(names(dots),names.legend)]
  if(type1=="l" | type1=="b") dots.leg <- c(dots.leg,list(lty=lty1))
  if(!(type1 %in% c("l","b"))) dots.leg <- dots[setdiff(names(dots),"lwd")]
    
  if(!is.null(legend)) do.call("legend",c(dots.leg,list(x=legend.loc, legend=legend, box.lty=box.lty1, inset=inset1,col=color)))
 }
 else
 {
  cat("Output can not be produced since the program stopped abnormally. \n")
  X <- NA
  Y <- NA
  lower <- NA
  upper <- NA
 }
 
 return(invisible(data.frame(time=X,Y,lower,upper))) 
}

plot.predict <- function(x,newdata,var.time,legend.loc="topright",na.action=1,legend,add=FALSE,confint=FALSE,...) UseMethod("plot.predict")
