print.Jointlcmm <-
function(x,...){
if (!inherits(x, "Jointlcmm")) stop("use only with \"Jointlcmm\" objects")

#cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")



cat("Joint latent class model for quantitative longitudinal outcome", "\n")
cat("     fitted by maximum likelihood method", "\n")

cl <- x$call
cat(" \n")
dput(cl)
cat(" \n")

cat("Statistical Model:", "\n")
cat(paste("     Dataset:", x$dataset),"\n")
cat(paste("     Number of subjects:", x$ns),"\n")
cat(paste("     Number of repeated measures:", length(x$pred[,1])),"\n")
cat(paste("     Number of events: ", x$N[7]," (",x$N[7]*100/x$ns,"%)"),"\n")
cat(paste("     Number of latents classes:", x$ng), "\n")
cat(paste("     Number of parameters:", length(x$best))," \n")


if (x$ng>1) {
if (x$hazard[[2]]=="Specific") cat("     Class-specific hazards and \n")
if (x$hazard[[2]]=="PH") cat("     Proportional hazards over latent classes  and \n")
if (x$hazard[[2]]=="Common") cat("  Common hazards over classes  and \n")
}

if (x$hazard[[1]]==2) {
cat("     Weibull baseline risk function \n")
}
if (x$hazard[[1]]==1)
{
cat("     Piecewise constant baseline risk function \n")
}
if (x$hazard[[1]]==3)
{
cat("     M-splines constant baseline risk function with nodes \n")
cat(     x$hazard[[3]]," \n")
}



cat(" \n")
cat("Iteration process:", "\n")

if(x$conv==1) cat("Convergence criteria satisfied")
if(x$conv==2) cat("Maximum number of iteration reached without convergence")
if(x$conv==4|x$conv==12){

cat("The program stopped abnormally. No results can be provided.\n")
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
cat(" \n")

}
}

