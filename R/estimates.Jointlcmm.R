
estimates.Jointlcmm <- function(x,cholesky=TRUE)
{
 if(missing(x)) stop("The argument x should be specified")
 if (!inherits(x, "Jointlcmm")) stop("use only with \"Jointlcmm\" objects")
 if(is.na(as.logical(cholesky))) stop("cholesky should be TRUE or FALSE")

 if(x$conv==1 | x$conv==2)
 {
  res <- x$best
  cholesky <- as.logical(cholesky)

  N <- x$specif[[1]]
  if(isTRUE(cholesky) & N[5]>0)
  {
   if(!isTRUE(x$call$idiag) | N[5]==1)
   {
    res[N[4]+1:N[5]] <- x$cholesky
   }
   else
   {
    nea <- sum(x$specif[[4]])
    res[N[4]+1:N[5]] <- x$cholesky[-setdiff(1:(nea*(nea+1)/2),1:nea*(1:nea+1)/2)]
   }

   names(res) <- sub("Varcov","cholesky",names(res))
  }
 }
 else
 {
  res <- NA
  cat("Output can not be produced since the program stopped abnormally. \n")
 }
 
 return(res)
}

estimates <- function(x,cholesky=TRUE) UseMethod("estimates")
