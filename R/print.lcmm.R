print.lcmm <-
function(x,...){
if (!inherits(x, "lcmm")) stop("use only with \"lcmm\" objects")

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
cat(paste("     Number of observations:", x$N[5]),"\n")
if(length(x$na.action))cat(paste("     Number of observations deleted:",length(x$na.action)),"\n")
cat(paste("     Number of latent classes:", x$ng), "\n")
cat(paste("     Number of parameters:", length(x$best))," \n")
if (x$linktype==0) {
ntrtot <- 1
cat("     Link function: linear"," \n")
}
if (x$linktype==1)
{
ntrtot <- 3
cat("     Link function: Standardised Beta CdF"," \n")
}
if (x$linktype==2) {
ntrtot <- length(x$linknodes)+1
cat("     Link function: Quadratic I-splines with nodes"," \n")
cat(     x$linknodes," \n")
}

cat(" \n")
cat("Iteration process:", "\n")

if(x$conv==1) cat("     Convergence criteria satisfied")
if(x$conv==2) cat("     Maximum number of iteration reached without convergence")
if(x$conv==4|x$conv==12) {
cat("     The program stopped abnormally. No results can be displayed.\n")
}else{

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
if (x$Ydiscrete==1){
cat(paste("     Discrete posterior log-likelihood:", round(x$discrete_loglik,2))," \n")
cat(paste("     Discrete AIC:", round(-2*x$discrete_loglik+2*length(x$best),2))," \n")
cat(" \n")

cat(paste("     Mean discrete AIC per subject:",round((-x$discrete_loglik+length(x$best))/as.double(x$ns),4))," \n")
cat(paste("     Mean UACV per subject:",round(x$UACV,4))," \n")
cat(paste("     Mean discrete LL per subject:",round(x$discrete_loglik/as.double(x$ns),4))," \n")
}
cat(" \n")
}}
