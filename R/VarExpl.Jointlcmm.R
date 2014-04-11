
VarExpl.Jointlcmm <- function(x,values)
{
 if(missing(x)) stop("The model should be specified")
 if (!inherits(x, "Jointlcmm")) stop("use only with \"Jointlcmm\" objects")
 if(missing(values)) values <- data.frame("intercept"=1)
 if (!inherits(values, "data.frame")) stop("values should be a data.frame object")
 if(any(is.na(values))) stop("values should not contain any missing values")


 if(x$conv==1 | x$conv==2)
 {
  res <- matrix(0,nrow=1,ncol=x$specif[[3]])

  N <- x$specif[[1]]
  Xnames <- x$Names[[2]]
  idea <- x$specif[[4]]
  
  names.random <- NULL
  name.cor <- NULL
  if(N[5]>0) names.random <- Xnames[which(idea==1)]
  if(length(N)>9 & N[10]>0) name.cor <- Xnames[which(x$specif[[10]]==1)]

  if(!is.null(names.random) | !is.null(name.cor))
  {
   names.values <- unique(c(names.random,name.cor))   #contient I(T^2))

   vars <- unique(c(all.vars(x$call$random),all.vars(x$call$cor)))
   if(!all(vars %in% colnames(values))) stop(paste(c("values should give a value for each of the following covariates: ","\n",vars,collapse=" ")))

   ### pour les facteurs

   #cas où une variable du dataset est un facteur
   olddata <- eval(x$call$data)
   for(v in setdiff(vars,"intercept"))
   {
    if(is.factor(olddata[,v]))
    {
     mod <- levels(olddata[,v])
     if (!(levels(as.factor(values[,v])) %in% mod)) stop(paste("invalid level in factor", v))
     values[,v] <- factor(values[,v], levels=mod)
    }
   }

   #cas où on a factor() dans l'appel
   call_random <- x$call$random
   z <- all.names(call_random)
   ind_factor <- which(z=="factor")
   if(length(ind_factor))
   {
    nom.factor <- z[ind_factor+1]
    for (v in nom.factor)
    {
     mod <- levels(as.factor(olddata[,v]))
     if (!all(levels(as.factor(values[,v])) %in% mod)) stop(paste("invalid level in factor", v))
     values[,v] <- factor(values[,v], levels=mod)
    }
   }
   call_random <- gsub("factor","",call_random)

   if(!is.null(name.cor)) values1 <- model.matrix(formula(paste(call_random[2],name.cor,sep="")),data=values)
   else values1 <- model.matrix(formula(call_random),data=values)

   if(colnames(values1)[1]=="(Intercept)") colnames(values1)[1] <- "intercept"

   if(nrow(values1)>1) warning("only the first line of values will be used")
   var.random <- values1[1,names.random]
   var.cor <- values1[1,name.cor]

   nea <- sum(idea==1)
   VarU <- matrix(0,nea,nea)
   if(nea==N[5])
   {
    diag(VarU) <- x$best[N[4]+1:N[5]]
   }
   else
   {
    VarU[lower.tri(VarU,diag=TRUE)] <- x$best[N[4]+1:N[5]]
    VarU <- t(VarU)
    VarU[lower.tri(VarU,diag=TRUE)] <- x$best[N[4]+1:N[5]]
   }

   numer <- t(var.random) %*% VarU %*% var.random
   if(x$specif[[3]]>0)
   {
    nw <- rep(1,x$specif[[3]])
    if(N[6]>0) nw <- c((x$best[N[4]+N[5]+1:N[6]])^2,1)
    numer <- numer * nw
   }

   Corr <- 0
   if(length(N)>9 & N[10]>0)
   {
    var.cor <- values[name.cor]
    if(N[10]==1)
    {
     Corr <- (x$best[N[4]+N[5]+N[6]])^2 * var.cor
    }
    if(N[10]==2)
    {
     Corr <- (x$best[N[4]+N[5]+N[6]])^2
    }
   }

   denom <- numer + Corr + (x$best[length(x$best)])^2

   res[1,] <- numer/denom *100
  }

  rownames(res) <- "%Var"
  colnames(res) <- paste("class",1:x$specif[[3]],sep="")
 }
 else
 {
  cat("Output can not be produced since the program stopped abnormally. \n")
  res <- NA
 }


 return(res)
}

VarExpl <- function(x,values) UseMethod("VarExpl")
