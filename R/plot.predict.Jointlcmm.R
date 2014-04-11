

plot.predict.Jointlcmm <- function(x,newdata,var.time,legend.loc="topright",na.action=1,legend,add=FALSE,confint=FALSE,...)
{
 if(missing(var.time)) stop("The argument var.time should be specified")
 if(missing(newdata)) stop("The argument newdata should be specified")
 if(missing(x)) stop("The argument x should be specified")
 if (!inherits(x, "Jointlcmm")) stop("use only with \"Jointlcmm\" objects")
 #--------------->
 if(!all(x$Names2[[1]] %in% c(colnames(newdata),"intercept"))) {
 stop(paste(c("newdata should at least include the following covariates: ","\n",x$Names2[[1]][-1]),collapse=" "))
 #cat("newdata should at least include the following covariates: ", "\n")
 #cat(x$Names2[[1]][-1], "\n")
 }
 #if (!all(x$Names2[[1]] %in% c(colnames(newdata),"intercept"))) stop("see above")
 #--------------->
 
 
 if (!inherits(newdata, "data.frame")) stop("newdata should be a data.frame object")
 if (!inherits(var.time, "character")) stop("the class of var.time should be character")
 if (!(var.time %in% colnames(newdata))) stop("The variable names var.time should be included in newdata")
 if(is.na(as.logical(add))) stop("add should be TRUE or FALSE")
   
 if(x$conv==1|x$conv==2)
 {
  newdata<-newdata[sort.list(newdata[,var.time]),] 
  
  #------------> changement Cecile 10/04/2012 ad
  if(x$Names2[[1]][1]!="intercept"){
  	newdata1 <- newdata[,x$Names2[[1]]]
  	colnames(newdata1) <- x$Names2[[1]]
  	newdata1 <- data.frame(newdata1)
  }else{
  	newdata1 <- cbind(intercept=rep(1,length=length(newdata[,1])),newdata[,x$Names2[[1]][-1]])
  	colnames(newdata1) <- x$Names2[[1]]
  	newdata1 <- data.frame(newdata1)
  }
  #------------>
  X1 <- NULL
  X2 <- NULL
  b1 <- NULL
  b2 <- NULL
  
  
  call_fixed <- x$call$fixed[3]
  if(is.null(x$call$random)) {call_random <- ~-1} else call_random <- x$call$random
  if(is.null(x$call$classmb)) {call_classmb <- ~-1} else call_classmb <- x$call$classmb
  if(is.null(x$call$mixture)) {call_mixture <- ~-1} else call_mixture <- x$call$mixture
  if(is.null(x$call$survival)) {call_survival <- ~-1} else call_survival <- x$call$survival[3]
  
  
  if(!(na.action%in%c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument") 
  
  if(na.action==1){
  	na.action=na.omit
  }else{
  	na.action=na.fail
  }
  
  ### pour les facteurs
  
  Xnames2 <- x$Names2[[1]]
  
   #cas où une variable du dataset est un facteur
   olddata <- eval(x$call$data)
    for(v in Xnames2[-1])
   {
    if (is.factor(olddata[,v]) & !(is.factor(newdata[,v])))
    {
     mod <- levels(olddata[,v])
     if (!(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
     newdata1[,v] <- factor(newdata1[,v], levels=mod)
    }
   }
   
   #cas où on a factor() dans l'appel
   z <- all.names(call_fixed)
   ind_factor <- which(z=="factor")
   if(length(ind_factor))
   {
    nom.factor <- z[ind_factor+1]  
    for (v in nom.factor)
    {
     mod <- levels(as.factor(olddata[,v]))
     if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
     newdata1[,v] <- factor(newdata1[,v], levels=mod)
    }
   }
   call_fixed <- gsub("factor","",call_fixed)
  
   z <- all.names(call_random)
   ind_factor <- which(z=="factor")
   if(length(ind_factor))
   {
    nom.factor <- z[ind_factor+1]
    for (v in nom.factor)
    {
     mod <- levels(as.factor(olddata[,v]))
     if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
     newdata1[,v] <- factor(newdata1[,v], levels=mod)
    }
   }
   call_random <- gsub("factor","",call_random)
         
   z <- all.names(call_classmb)
   ind_factor <- which(z=="factor")
   if(length(ind_factor))
   {
    nom.factor <- z[ind_factor+1]
    for (v in nom.factor)
    {
     mod <- levels(as.factor(olddata[,v]))
     if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
     newdata1[,v] <- factor(newdata1[,v], levels=mod)
    }
   }
   call_classmb <- gsub("factor","",call_classmb)
          
   z <- all.names(call_mixture)
   ind_factor <- which(z=="factor")
   if(length(ind_factor))
   {
    nom.factor <- z[ind_factor+1]
    for (v in nom.factor)
    {
     mod <- levels(as.factor(olddata[,v]))
     if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
     newdata1[,v] <- factor(newdata1[,v], levels=mod)
    }
   }
   call_mixture <- gsub("factor","",call_mixture) 
  
   z <- all.names(call_survival)
   ind_factor <- which(z=="factor")
   if(length(ind_factor))
   {
    nom.factor <- z[ind_factor+1]  
    for (v in nom.factor)
    {
     mod <- levels(as.factor(olddata[,v]))
     if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
     newdata1[,v] <- factor(newdata1[,v], levels=mod)
    }
   }
   call_survival <- gsub("factor","",call_survival)
   call_survival <- gsub("mixture","",call_survival)
   
     
       
  ### Traitement des donnees manquantes
  # fixed
  #mcall <- x$call[c(1,match(c("data"),names(x$call),0))]
  mcall <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
  mcall$na.action <- na.action
  mcall$data <- newdata1
  
  
  # fixed
  m <- mcall
  m$formula <- formula(paste("~",call_fixed,sep=""))
  m[[1]] <- as.name("model.frame")	
  m <- eval(m, sys.parent()) 
  na.fixed <- attr(m,"na.action")
  
  # mixture
  if(!is.null(x$call$mixture)){
  	m <- mcall
  	m$formula <- formula(paste("~",call_mixture,sep=""))
  	m[[1]] <- as.name("model.frame")	
  	m <- eval(m, sys.parent()) 
  	na.mixture <- attr(m,"na.action")
  }else{
  	na.mixture <- NULL
  }
  
  # random
  if(!is.null(x$call$random)){
  	m <- mcall
  	m$formula <- formula(paste("~",call_random,sep=""))
  	m[[1]] <- as.name("model.frame")	
  	m <- eval(m, sys.parent()) 
   	na.random <- attr(m,"na.action")
  }else{
  	na.random <- NULL
  }
  
  # classmb
  if(!is.null(x$call$classmb)){ 
  	m <- mcall	
  	m$formula <- formula(paste("~",call_classmb,sep=""))
  	m[[1]] <- as.name("model.frame")	
  	m <- eval(m, sys.parent()) 
   	na.classmb <- attr(m,"na.action")
  }else{
  	na.classmb <- NULL
  }
  
  #survival
  if(!is.null(x$call$survival))
  {
   m <- mcall
   m$formula <- formula(paste("~",call_survival,sep=""))
   m[[1]] <- as.name("model.frame")	
   m <- eval(m, sys.parent()) 
   na.survival <- attr(m,"na.action")
  }
  else {na.survival <- NULL}
  
  ########### For survival
  #if(!is.null(x$call$survival)){ 
  #	res.evt <- x$call$survival
  #	class(res.evt) <- "formula"
  #	tmp.res <- res.evt
  #
  #	res.evt <- terms(res.evt,"mixture") 
  #	inddep.surv <- attr(res.evt, "term.labels")
  #	ind.mixture <- untangle.specials(res.evt, "mixture", 1)
  #	inddepvar.Mixt  <- gsub("\\)","",gsub("mixture\\(","",ind.mixture$vars))
  #	inddepvar.noMixt <- inddep.surv[!(inddep.surv %in% ind.mixture$vars)]
  #	tmp.su <- c(inddepvar.noMixt,inddepvar.Mixt)
  #	names.survival <- "~" 
  #	for(i in 1:length(tmp.su)){
  #	      if(i==1){
  #		      names.survival <- paste(names.survival,tmp.su[i],sep="")
  #	      }else{
  #		      names.survival <- paste(names.survival,tmp.su[i],sep="+")
  #	      }
  #	}
  
  ##7/05/2012
  ##7/05/2012
  #	m <- mcall
  #	m$formula <- formula(names.survival)
  #	m[[1]] <- as.name("model.frame")
  #	m <- eval(m, sys.parent())
  #	na.survival <- attr(m,"na.action")
  
  ## Table sans donnees manquante: newdata
  	na.action <- unique(c(na.fixed,na.mixture,na.random,na.classmb,na.survival))
  	if(!is.null(na.action)){
  		newdata1 <- newdata1[-na.action,]
  	}
  ##7/05/2012
  #	X_survival <- model.matrix(formula(names.survival),data=newdata1)
  #	if(colnames(X_survival)[1]=="(Intercept)"){
  #		colnames(X_survival)[1] <- "intercept"
  #	}
  #	id.X_survival <- 1
  #}else{
  #	id.X_survival <- 0
  ##7/05/2012
  #	na.action <- unique(c(na.fixed,na.mixture,na.random,na.classmb))
  #	if(!is.null(na.action)){
  #		newdata1 <- data[-na.action,]
  #	}
  #}
  
  
  
  ## Construction de nouvelles var explicatives sur la nouvelle table
  ## fixed
  	
  	X_fixed <- model.matrix(formula(paste("~",call_fixed,sep="")),data=newdata1)
  	if(colnames(X_fixed)[1]=="(Intercept)"){
  		colnames(X_fixed)[1] <- "intercept"
  		int.fixed <- 1
  	}	
  ## mixture
  	if(!is.null(x$call$mixture)){
  		X_mixture <- model.matrix(formula(paste("~",call_mixture,sep="")),data=newdata1)	
  		if(colnames(X_mixture)[1]=="(Intercept)"){
  			colnames(X_mixture)[1] <- "intercept"
  			int.mixture <- 1
  		}
  		id.X_mixture <- 1
  	}else{
  		id.X_mixture <- 0
  	}	
  ## random
  	if(!is.null(x$call$random)){
  		X_random <- model.matrix(formula(paste("~",call_random,sep="")),data=newdata1)	
  		if(colnames(X_random)[1]=="(Intercept)"){
  			colnames(X_random)[1] <- "intercept"
  			int.random <- 1
  		}
  		id.X_random <- 1
  	}else{
  		id.X_random <- 0
  	}	
  ## classmb
  	if(!is.null(x$call$classmb)){ 
  		X_classmb <- model.matrix(formula(paste("~",call_classmb,sep="")),data=newdata1)
  		colnames(X_classmb)[1] <- "intercept"
  		id.X_classmb <- 1
  	}else{
  		id.X_classmb <- 0
  	}	
  
  ## survival
  	if(!is.null(x$call$survival)){ 
  		X_survival <- model.matrix(formula(paste("~",call_survival,sep="")),data=newdata1)
  		colnames(X_survival)[1] <- "intercept"
  		id.X_survival <- 1
  	}else{
  		id.X_survival <- 0
  	}	
  	
  	
  ## Construction des var expli
  newdata1 <- X_fixed
  
  if(id.X_mixture == 1){
  	for(i in 1:length(colnames(X_mixture))){
  		if((colnames(X_mixture)[i] %in% colnames(newdata1))==F){
  			newdata1 <- cbind(newdata1,X_mixture[,i])
  			
  		}
  	}
  }
  if(id.X_random == 1){
  	for(i in 1:length(colnames(X_random))){
  		if((colnames(X_random)[i] %in% colnames(newdata1))==F){
  			newdata1 <- cbind(newdata1,X_random[,i])
  		}	 
  	}
  }
  if(id.X_classmb == 1){
  	for(i in 1:length(colnames(X_classmb))){
  		if((colnames(X_classmb)[i] %in% colnames(newdata1))==F){
  			newdata1 <- cbind(newdata1,X_classmb[,i],deparse.level=0)	 
  		}	
  	}
  }
  if(id.X_survival == 1){
  	for(i in 1:length(colnames(X_survival))){
  		if((colnames(X_survival)[i] %in% colnames(newdata1))==F){
  			newdata1 <- cbind(newdata1,X_survival[,i])
  			
  		}
  	}
  }
  
  #colnames(newdata1) <- namesNew
  ### end add 11/04/2012
  
  placeV <- list() #places pour les variances
  placeV$commun <- NA
  for(i in 1:x$specif[[3]])
  {
   placeV[paste("class",i,sep="")] <- NA
  }  
  
  kk<-0
  for(k in 1:length(x$specif[[6]]))
  {
   if(x$specif[[6]][k]==1)
   {
    X1 <- cbind(X1,newdata1[,k])
    place <- x$specif[[1]][1]+x$specif[[1]][2]+x$specif[[1]][3]+kk
    b1 <- c(b1,x$best[place+1])
    placeV$commun <- c(placeV$commun,place+1)  
    kk <- kk+1
   }
  
   if(x$specif[[6]][k]==2)
   {
    X2 <- cbind(X2,newdata1[,k])
    place1 <- x$specif[[1]][1]+x$specif[[1]][2]+x$specif[[1]][3]+kk+1
    place2 <- x$specif[[1]][1]+x$specif[[1]][2]+x$specif[[1]][3]+kk+x$specif[[3]]
    b2 <- rbind(b2,x$best[place1:place2])
    for(i in 1:x$specif[[3]])
    {
     placeV[[paste("class",i,sep="")]] <- c(placeV[[paste("class",i,sep="")]],x$specif[[1]][1]+x$specif[[1]][2]+x$specif[[1]][3]+kk+i)
    }
    kk <- kk+x$specif[[3]]
   }
  }
  
  Y<-matrix(0,length(newdata1[,1]),x$specif[[3]])
  colnames(Y) <- paste("class",1:x$specif[[3]],sep="") 
  for(g in 1:x$specif[[3]])
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
  #ylim1 <- ylim
  #if (is.null(ylim)){
  #ylim1 <- c(min(Y),max(Y))
  #}
  
  #extraction de Var(beta)
  nbeta <- x$specif[[1]][4]-x$specif[[1]][1]-x$specif[[1]][2]-x$specif[[1]][3]
  Vbeta <- matrix(0,nbeta,nbeta)
  npm <- length(x$best)
  indice <- 1:npm * (1:npm+1) /2
  indtmp <- indice[(x$specif[[1]][4]-nbeta+1):x$specif[[1]][4]]
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
  colnames(lower) <- paste("lower.class",1:x$specif[[3]],sep="")
  colnames(upper) <- paste("upper.class",1:x$specif[[3]],sep="") 
   
  if(x$specif[[3]]==1)
  { 
   varpred <- apply(X1,1,function(x) matrix(x,nrow=1) %*% Vbeta %*% matrix(x,ncol=1))  #controler si l'ordre est le meme entre X1 et Vbeta
   lower[,1] <- Y[,1] -1.96 * sqrt(varpred)
   upper[,1] <- Y[,1] +1.96 * sqrt(varpred)
  }
  else
  {
   for(g in 1:x$specif[[3]])
   {
    ind <- na.omit(c(placeV[["commun"]],placeV[[paste("class",g,sep="")]]))
    X12 <- cbind(X1,X2)
    X12 <- X12[,order(ind)]
    
    varclass <- Vbeta[sort(ind)-sum(x$specif[[1]][1:3]),sort(ind)-sum(x$specif[[1]][1:3])]
    varpred <- diag(X12 %*% varclass %*% t(X12))
    
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
    else  color <- 1:x$specif[[3]]
  
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
    
    if(missing(legend)) legend <- paste("class",1:x$specif[[3]],sep="")
  
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
  
  #for (i in 1:x$specif[[3]]){
  #plot(Y[,i]~X,col=i,type="l",ylim=ylim1,xlim=c(min(X),max(X)),xlab=var.time,pch=2,ylab="outcome",main=title1,lty=i,...)
  #par(new=TRUE)
  #}
  #legend(x=legend.loc,legend=colnames(Y),col=1:x$specif[[3]],lty=1:x$specif[[3]],inset=.02,bty="n")
  #par(new=FALSE)
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
