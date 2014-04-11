dynpred <- function(model,newdata,landmark,horizon,var.time,
           fun.time=identity,na.action=1,draws=FALSE,ndraws=2000)    
{                                                       
 if(missing(model)) stop("The argument model must be specified")
 if(class(model)!="Jointlcmm") stop("The argument model must be a 'Jointlcmm' object")
 if(missing(newdata)) stop("The argument newdata should be specified")
 if (!inherits(newdata, "data.frame")) stop("newdata should be a data.frame object")
 if(missing(landmark)) stop("Please specify at least one landmark time")
 if(missing(horizon)) stop("Please specify at least one horizon time")
 if(any(horizon<=0)) stop("horizon must be positive times")
 if(missing(var.time)) stop("The argument 'var.time' should be specified")
 if(!is.character(var.time)) stop("'var.time' should be a character")
 if(!(var.time %in% colnames(newdata))) stop("'var.time' should be a variable included in 'newdata'") 
 if(!is.function(fun.time)) stop("'fun.time' should be a function")

#if (!all(x$Xnames2 %in% colnames(newdata))) stop(paste(c("newdata should at least include the following covariates: ","\n",x$Xnames2),collapse=" "))
 
 if(model$conv==1|model$conv==2)
 {
  if(model$conv==2 & draws==TRUE)
  {
   cat("No confidence interval will be provided since the program did not converge properly \n")
   draws <- FALSE
  }
  
  nbland <- length(landmark)
  nbhoriz <- length(horizon) 
  
  idprob <- model$specif[[5]]
  idea <- model$specif[[4]]
  idg <- model$specif[[6]]
  idcor <- model$specif[[10]]
  idxevt <- model$specif[[7]]
  idiag <- model$specif[[8]]
  nv <- length(idprob)
  nwg <- model$specif[[1]][6]
  ng <- model$specif[[3]]
  ncor <- model$specif[[1]][10]
  nz <- model$specif[[1]][8]
  zi <- model$hazard[[3]]
  typrisq <-  model$hazard[[1]]
  risqcom <- switch(model$hazard[[2]],"Specific"=0,"PH"=2,"Common"=1)
  logspecif <- model$specif[[9]]
  nvdepsurv <- model$specif[[1]][9] # 1 si une variable TimeDepVar , 0 sinon
  nvarxevt <- nvdepsurv + sum(idxevt==1) + ng*sum(idxevt==2)
  best <- model$best
  npm <- length(best) 
  
   #mettre cholesky a la place de varcov des effets aleatoires
   if(model$specif[[1]][5]>0)
   {   
    best[model$specif[[1]][4]+1:model$specif[[1]][5]] <- na.omit(model$cholesky)
   }
   
  call_fixed <- model$call$fixed[3]
  if(is.null(model$call$random)) {call_random <- ~-1} else call_random <- model$call$random
  if(is.null(model$call$classmb)) {call_classmb <- ~-1} else call_classmb <- model$call$classmb
  if(is.null(model$call$mixture)) {call_mixture <- ~-1} else call_mixture <- model$call$mixture
  if(is.null(model$call$survival)) {call_survival <- ~-1} else call_survival <- model$call$survival[3]


  if(!(na.action%in%c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument")

  if(na.action==1)
  {
  	na.action=na.omit
  }
  else
  {
  	na.action=na.fail
  }

  modelNames <- c(model$Names2[[1]],   # variables explicatives (MM et survie)
                  model$Names[[1]],    # nom de l'outcome du MM
                  model$Names[[5]],    # identifiant sujets
                  model$Names[[4]],    # nom prior
                  model$Names[[7]])    # nom TimeDepVar
                  
  if(length(model$Names[[3]])>2)
  {
   modelNames <- c(modelNames,model$Names[[3]][1]) #Tentry
  }                

 # if(model$Names2[[1]][1]!="intercept")
 # {
 # 	newdata1 <- newdata[,modelNames,drop=FALSE]
 # 	colnames(newdata1) <- modelNames
  #	newdata1 <- data.frame(newdata1)
  #}
  #else
  #{
 # 	newdata1 <- cbind(intercept=rep(1,length=length(newdata[,1])),newdata[,modelNames[-1],drop=FALSE])
 # 	colnames(newdata1) <- modelNames
 # 	newdata1 <- data.frame(newdata1)
 # }
  if(model$Names2[[1]][1]=="intercept")
  {
   newdata1 <- cbind(intercept=rep(1,length=length(newdata[,1])),newdata[,setdiff(colnames(newdata),model$Names2[[1]][1]),drop=FALSE])
   colnames(newdata1) <- c("intercept",setdiff(colnames(newdata),model$Names2[[1]][1]))
   newdata1 <- data.frame(newdata1)
  } 
  
  ### faire ici verif des variables dans newdata1
   #mettre peut-etre des valeurs par defaut dans prior et TimedepVar ?
  if(!all(modelNames %in% colnames(newdata1))) stop(paste(c("newdata should at least include the following covariates: ","\n",modelNames),collapse=" "))

  ### ordonner selon numero et temps
  newdata1 <- newdata1[order(newdata1[,model$Names[[5]]],newdata1[,var.time]),,drop=FALSE]  
    #est ce que fun.time peut modifier l'ordre????

  ### pour les facteurs

  Xnames2 <- model$Names2[[1]]

   #cas o? une variable du dataset est un facteur
   olddata <- eval(model$call$data)
   for(v in Xnames2[-1])
   {
    if (is.factor(olddata[,v]) & !(is.factor(newdata[,v])))
    {
     mod <- levels(olddata[,v])
     if (!(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
     newdata1[,v] <- factor(newdata1[,v], levels=mod)
    }
   }

   #cas o? on a factor() dans l'appel
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
  if(!is.null(model$call$mixture))
  {
   m <- mcall
   m$formula <- formula(paste("~",call_mixture,sep=""))
   m[[1]] <- as.name("model.frame")
   m <- eval(m, sys.parent())
   na.mixture <- attr(m,"na.action")
  }
  else
  {
   na.mixture <- NULL
  }

  # random
  if(!is.null(model$call$random))
  {
   m <- mcall
   m$formula <- formula(paste("~",call_random,sep=""))
   m[[1]] <- as.name("model.frame")
   m <- eval(m, sys.parent())
   na.random <- attr(m,"na.action")
  }
  else
  {
   na.random <- NULL
  }

  # classmb
  if(!is.null(model$call$classmb))
  {
   m <- mcall
   m$formula <- formula(paste("~",call_classmb,sep=""))
   m[[1]] <- as.name("model.frame")
   m <- eval(m, sys.parent())
   na.classmb <- attr(m,"na.action")
  }
  else
  {
   na.classmb <- NULL
  }

  #survival
  if(!is.null(model$call$survival))
  {
   m <- mcall
   m$formula <- formula(paste("~",call_survival,sep=""))
   m[[1]] <- as.name("model.frame")
   m <- eval(m, sys.parent())
   na.survival <- attr(m,"na.action")
  }
  else
  {
   na.survival <- NULL
  }
  
  #cor
  na.cor <- NULL
  if(length(model$specif[[1]])>9)
  {
   if(model$specif[[1]][10]>0)
   {
    z <- which(model$specif[[10]]==1)
    var.cor <- newdata1[,model$Names[[2]][z]]
    na.cor <- which(is.na(var.cor))
   }
  }
  
  #var.time
  na.time <- NULL
  if(!(var.time %in% model$Names2[[1]]))
  {
   na.time <- which(is.na(newdata1[,var.time]))
  }
  
  # toutes les variables pas var expl
  dtmp <- newdata1[,setdiff(modelNames,model$Names2[[1]]),drop=FALSE]

  m <- model.frame(formula=formula(paste(model$Names[[1]],"~",paste(setdiff(modelNames,model$Names2[[1]]),collapse="+"),sep="")),
  data=dtmp,na.action=na.action)
  na.other <- attr(m,"na.action")
  
  # fun.time 
  if(length(match.call()$fun.time))
  {   
   if(!isTRUE(as.character(match.call()[[which(names(match.call())=="fun.time")]])=="identity"))
   {
    z <- match.call()$fun.time 
   
    if(isTRUE(as.character(z) %in% ls(.GlobalEnv)))
    {
     bodyfun <- as.expression(body(get(as.character(z))))
    }
    else
    {
     bodyfun <- as.expression(z[[3]])
    }
    vars <- intersect(colnames(newdata1),all.names(bodyfun))
    getvars <- paste("getElement(newdata1,'",vars,"')",sep="")
 
    if(length(vars))
    {
     for(i in 1:length(vars))
     {
      bodyfun <- sub(vars[i],getvars[i],bodyfun)
      bodyfun <- sub(as.character(match.call()$newdata),"",bodyfun)
      tmp <- strsplit(bodyfun,split="")[[1]]
      dollard <- which(tmp=="$")
      if(length(dollard))
      {
       tmp <- tmp[-dollard]
      }
      bodyfun <- paste(tmp,collapse="")
     }
    }

    headfun <- paste("function(",names(formals(fun.time)),")",collapse="")
    ff1 <- paste(headfun,bodyfun)
    ff2 <- parse(text=ff1)
    ff3 <- eval(ff2)

    timemes <- do.call("ff3",list(newdata1[,var.time]))
   }
   else
   {
    timemes <- newdata1[,var.time]
   } 
  }
  else
  {
   timemes <- newdata1[,var.time]
  }
  na.fun <- which(is.na(timemes))
  
  ## Table sans donnees manquante: newdata1
  na.action <- unique(c(na.other,na.fixed,na.mixture,na.random,na.classmb,na.survival,na.cor,na.time,na.fun))
  if(length(na.action))
  {
   newdata1 <- newdata1[-na.action,]
  }
 	
  ## si Tevent et Devent dans newdata
  Tevent <- NULL
  Devent <- NULL
  
  if(length(model$Names[[3]])==2)  # ordre : Tevent,Event 
  {
   if(all(model$Names[[3]][1:2] %in% colnames(newdata1)))
   {
    Tevent <- newdata1[,model$Names[[3]][1]]
    Devent <- newdata1[,model$Names[[3]][2]]
    
    if(length(na.action))
    {
      Tevent <- Tevent[-na.action]
      Devent <- Devent[-na.action]
    }
    
    indNA <- which(is.na(Tevent) | is.na(Devent))
    if(length(indNA)) 
    {
     Tevent[indNA] <- 0
     Devent[indNA] <- 0
    } 
   }
  }

  if(length(model$Names[[3]])==3)  # ordre : Tentry,Tevent,Event 
  {
   if(all(model$Names[[3]][2:3] %in% colnames(newdata1)))
   {
    Tevent <- newdata1[,model$Names[[3]][2]]
    Devent <- newdata1[,model$Names[[3]][3]]
  
    if(!is.null(na.action))
    {
     Tevent <- Tevent[-na.action]
     Devent <- Devent[-na.action]
    }  
  
    indNA <- which(is.na(Tevent) | is.na(Devent))
    if(length(indNA)) 
    {
     Tevent[indNA] <- 0
     Devent[indNA] <- 0
    }  
   }
  }
              
  ## nb de sujets 
  id.subject <- newdata1[,model$Names[[5]]]
  ns <- length(unique(id.subject))
    
  ## vecteur Y
  Y <- newdata1[,model$Names[[1]]]
  
  ## nb d'observations
  nobs <- length(Y)
  
  ## nb de mesures par sujet
  nmes <- as.vector(table(id.subject))
  
  ## prior
  prior <- rep(0,ns)
  if(length(model$Names[[4]]))
  {
   prior <- newdata1[cumsum(nmes),model$Names[[4]]]
  }
  
  ## age entree si type = counting
  tsurv0 <- rep(0,ns)
  idtrunc <- 0
  if(length(model$Names[[3]])>2)
  {
   idtrunc <- 1
   tsurv0 <- newdata1[cumsum(nmes),model$Names[[3]][1]] 
  }
  
  ## TimeDepVar
  tsurvint <- rep(max(landmark)+1,ns)
  #ind_survint <- rep(0,ns)
  if(length(model$Names[[7]]))
  {
   tsurvint <- newdata1[cumsum(nmes),model$Names[[7]]]
   #ind_survint[tsurvint<times[1]] <- 1
  }
  #attention : dans Jointlcmm, si tsurvint=NA => tsurvint=tevent et ind_survint=0
  #               -> ne pas prendre tsurvint dans na.action?
  
  ## reduire Tevent et Devent au nb de sujets
  if(!is.null(Tevent))
  {
   Tevent <- Tevent[cumsum(nmes)]
   Devent <- Devent[cumsum(nmes)]   
  } 


  ## Construction de nouvelles var explicatives sur la nouvelle table :
  
  ## fixed

  X_fixed <- model.matrix(formula(paste("~",call_fixed,sep="")),data=newdata1)
  if(colnames(X_fixed)[1]=="(Intercept)")
  {
   colnames(X_fixed)[1] <- "intercept"
   int.fixed <- 1
  }
  
  ## mixture
  if(!is.null(model$call$mixture))
  {
   X_mixture <- model.matrix(formula(paste("~",call_mixture,sep="")),data=newdata1)
   if(colnames(X_mixture)[1]=="(Intercept)")
   {
    colnames(X_mixture)[1] <- "intercept"
    int.mixture <- 1
   }
   id.X_mixture <- 1
  }
  else
  {
   id.X_mixture <- 0
  }
  	
  ## random
  if(!is.null(model$call$random))
  {
   X_random <- model.matrix(formula(paste("~",call_random,sep="")),data=newdata1)
   if(colnames(X_random)[1]=="(Intercept)")
   {
    colnames(X_random)[1] <- "intercept"
    int.random <- 1
   }
   id.X_random <- 1
  }
  else
  {
   id.X_random <- 0
  }
  
  ## classmb
  if(!is.null(model$call$classmb))
  {
   X_classmb <- model.matrix(formula(paste("~",call_classmb,sep="")),data=newdata1)
   colnames(X_classmb)[1] <- "intercept"
   id.X_classmb <- 1
  }
  else
  {
   id.X_classmb <- 0
  }

  ## survival
  if(!is.null(model$call$survival))
  {
   X_survival <- model.matrix(formula(paste("~",call_survival,sep="")),data=newdata1)
   colnames(X_survival)[1] <- "intercept"
   id.X_survival <- 1
  }
  else
  {
   id.X_survival <- 0
  }

  ##cor
  if(length(model$specif[[1]])>9)
  {
   if(model$specif[[1]][10]>0)  #on reprend la variable de temps de cor
   {
    z <- which(model$specif[[10]]==1)
    var.cor <- newdata1[,model$Names[[2]][z]]
   }
  }
  
  ## var.time
  if(length(na.action)) timemes <- timemes[-na.action]
  #  temps deja dans la bonne echelle
       

  ## Construction de newdata1 dans le bon ordre
  X <- X_fixed

  if(id.X_mixture == 1)
  {
   for(i in 1:length(colnames(X_mixture)))
   {
    if((colnames(X_mixture)[i] %in% colnames(X))==FALSE)
    {
     X <- cbind(X,X_mixture[,i])
    }
   }
  }
  if(id.X_random == 1)
  {
   for(i in 1:length(colnames(X_random)))
   {
    if((colnames(X_random)[i] %in% colnames(X))==FALSE)
    {
     X <- cbind(X,X_random[,i])
    }
   }
  }
  if(id.X_classmb == 1)
  {
   for(i in 1:length(colnames(X_classmb)))
   {
    if((colnames(X_classmb)[i] %in% colnames(X))==FALSE)
    {
     X <- cbind(X,X_classmb[,i],deparse.level=0)
    }
   }
  }
  if(id.X_survival == 1)
  {
   for(i in 1:length(colnames(X_survival)))
   {
    if((colnames(X_survival)[i] %in% colnames(X))==FALSE)
    {
     X <- cbind(X,X_survival[,i])
    }
   }
  }

  if(length(model$specif[[1]])>9)
  {
   if(model$specif[[1]][10]>0)
   {
    if(model$specif[[4]][z]==0 & model$specif[[5]][z]==0 & model$specif[[6]][z]==0 & model$specif[[7]][z]==0)
    {
     X <- cbind(X,var.cor)
    }
   }
  }


  
  if(!isTRUE(draws))
  {
   proba <- rep(0,ns*nbland*nbhoriz)
   
   out <- .Fortran("proba_evt_dyn",as.double(Y),as.double(X),as.integer(ns),as.integer(nmes),as.integer(nobs),
   as.integer(prior),as.integer(ng),as.integer(nv),as.integer(idiag),as.integer(nwg),as.integer(ncor),
   as.integer(logspecif),as.integer(idxevt),as.double(zi),as.integer(idea),as.integer(idg),as.integer(idprob),
   as.integer(idcor),as.integer(risqcom),as.integer(nvdepsurv),as.integer(nvarxevt),as.integer(typrisq),
   as.integer(nz),as.double(tsurv0),as.double(tsurvint),as.integer(idtrunc),
   as.double(best),as.integer(npm),as.double(timemes),as.double(landmark),as.double(horizon),
   as.integer(nbland),as.integer(nbhoriz),proba=as.double(proba),PACKAGE="lcmm")
 
 
   proba <- out$proba
      
   #mettre NA si proba=-5 (ie nmes=0)
   ind5 <- which(proba==-5)
   if(length(ind5)) proba[ind5] <- NA

   res <- rep(sort(unique(id.subject)),each=nbland*nbhoriz)
   res <- cbind(res,rep(rep(landmark,each=nbhoriz),ns))
   res <- cbind(res,rep(horizon,ns*nbland))
   res <- cbind(res,proba)
   colnames(res) <- c(model$Names[[5]],"landmark","horizon","pred")
    
   # mettre NA si Devent=1 et Tevent<s
   if(!is.null(Tevent))
   {
    matevt <- cbind(rep(Tevent,each=nbland*nbhoriz),rep(Devent,each=nbland*nbhoriz))
    indevt <- which( (matevt[,2]==1) & (matevt[,1]<res[,"landmark"]) )
    if(length(indevt)) res[indevt,"pred"] <- NA
   }                                                                              	
  }
  else   # ie avec draws
  {
   Mat <- matrix(0,ncol=length(best),nrow=length(best))
   Mat[upper.tri(Mat,diag=TRUE)]<- model$V
   Chol <- chol(Mat)
   Chol <- t(Chol) 	
   
   doOneDraw <- function()
   {
    bdraw <- rnorm(npm)
    bdraw <- best + Chol %*% bdraw
    
    proba <- rep(0,ns*nbland*nbhoriz)
    
    out <- .Fortran("proba_evt_dyn",as.double(Y),as.double(X),as.integer(ns),as.integer(nmes),as.integer(nobs),
    as.integer(prior),as.integer(ng),as.integer(nv),as.integer(idiag),as.integer(nwg),as.integer(ncor),
    as.integer(logspecif),as.integer(idxevt),as.double(zi),as.integer(idea),as.integer(idg),as.integer(idprob),
    as.integer(idcor),as.integer(risqcom),as.integer(nvdepsurv),as.integer(nvarxevt),as.integer(typrisq),
    as.integer(nz),as.double(tsurv0),as.double(tsurvint),as.integer(idtrunc),
    as.double(bdraw),as.integer(npm),as.double(timemes),as.double(landmark),as.double(horizon),
    as.integer(nbland),as.integer(nbhoriz),proba=as.double(proba),PACKAGE="lcmm")
      
    proba <- out$proba
   
    pb <- 0
    if(any(!is.finite(proba))) pb <- 1
    
    #mettre NA si proba=-5 (ie nmes=0)
    ind5 <- which(proba==-5)
    if(length(ind5)) proba[ind5] <- NA    

    return(c(pb,proba)) 
   }
  
   ndraws <- as.integer(ndraws)
 
   resdraws <- replicate(ndraws,doOneDraw())

   probs <- resdraws[-1,,drop=FALSE]
   pb <- sum(resdraws[1,])
  
   if(pb>0) warning("Infinite probabilities have been found. Confidence intervals are based on ", ndraws-pb ," simulations.")
  
   med <- apply(probs,1,median,na.rm=TRUE)
   qmin <- apply(probs,1,quantile,prob=0.025,na.rm=TRUE)
   qmax <- apply(probs,1,quantile,prob=0.975,na.rm=TRUE)
  
 
   res <- rep(sort(unique(id.subject)),each=nbland*nbhoriz)
   res <- cbind(res,rep(rep(landmark,each=nbhoriz),ns))
   res <- cbind(res,rep(horizon,ns*nbland))
   res <- cbind(res,med,qmin,qmax)
     
   colnames(res) <- c(model$Names[[5]],"landmark","horizon","pred","pred_2.5","pred_97.5")
  
   # mettre NA si Devent=1 et Tevent<s
   if(!is.null(Tevent))
   {
    matevt <- cbind(rep(Tevent,each=nbland*nbhoriz),rep(Devent,each=nbland*nbhoriz))
    indevt <- which( (matevt[,2]==1) & (matevt[,1]<res[,"landmark"]) )
    if(length(indevt)) res[indevt,c("pred","pred_2.5","pred_97.5")] <- NA
   }
  }
 }
 else  # ie conv != 1 ou 2
 {
  cat("Output can not be produced since the program stopped abnormally. \n")
  res <- NA
  id.subject <- NA
  nmes <- NA
  Y <- NA
  timemes <- NA 
 }
  
 prmatrix(res)  
 
 res.list <- list(pred=res,newdata=data.frame(id=rep(unique(id.subject),nmes)[which(timemes<max(landmark))],y=Y[which(timemes<max(landmark))],time=timemes[which(timemes<max(landmark))]))
 class(res.list) <- "dynpred"
 return(invisible(res.list))
}




