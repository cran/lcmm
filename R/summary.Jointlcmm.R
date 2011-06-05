summary.Jointlcmm <-
function(object,...){
x <- object
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
if(x$conv==4|x$conv==12) {
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
cat(paste("     Score test statistic for CI assumption: ", round(x$CIstat,3)," (p-value=",round((1-pchisq(x$CIstat,sum(x$idea0))),4),")" ,sep=""))
cat(" \n")
cat(" \n")
cat("Maximum Likelihood Estimates:", "\n")
cat(" \n")

NPROB <- x$N[1]
nrisqtot <-x$N[2]
nvarxevt <-x$N[3]
NEF   <- x$N[4]
NVC   <- x$N[5]
NW    <- x$N[6]
NPM   <- length(x$best)


se <- rep(1,NPM)
if ((all.equal(x$conv,1)==T)==T){
#recuperation des indices de V
id <- 1:NPM
indice <- rep(id*(id+1)/2)
se <-sqrt(x$V[indice])
se[(NEF+1):(NEF+NVC)]<-1

wald <- x$best/se
pwald <- 1-pchisq(wald**2,1)
coef <- x$best
}else{
se <- NA
wald <- NA
pwald <- NA
coef <- x$best
}





if(NPROB>0){
cat(" *** Fixed effects in the class-membership model:\n" )

tmp <- cbind(coef[1:NPROB],se[1:NPROB],wald[1:NPROB],pwald[1:NPROB])
dimnames(tmp) <- list(names(coef)[1:NPROB], c("coef", "Se", "Wald", "p-value"))
cat("\n")
prmatrix(tmp)
cat("\n")

}

if((nrisqtot+nvarxevt)>0){
cat(" *** Parameters in the proportional hazard model:\n" )

tmp <- cbind(coef[(NPROB+1):(NPROB+nrisqtot+nvarxevt)],se[(NPROB+1):(NPROB+nrisqtot+nvarxevt)],wald[(NPROB+1):(NPROB+nrisqtot+nvarxevt)],pwald[(NPROB+1):(NPROB+nrisqtot+nvarxevt)])
dimnames(tmp) <- list(names(coef)[(NPROB+1):(NPROB+nrisqtot+nvarxevt)], c("coef", "Se", "Wald", "p-value"))
cat("\n")
prmatrix(tmp)
cat("\n")

}


cat(" *** Fixed effects in the longitudinal model:\n" )

tmp <- cbind(coef[(NPROB+nrisqtot+nvarxevt+1):(NEF)],se[(NPROB+nrisqtot+nvarxevt+1):(NEF)],wald[(NPROB+nrisqtot+nvarxevt+1):(NEF)],pwald[(NPROB+nrisqtot+nvarxevt+1):(NEF)])
dimnames(tmp) <- list(names(coef)[(NPROB+nrisqtot+nvarxevt+1):(NEF)], c("coef", "Se", "Wald", "p-value"))
cat("\n")
prmatrix(tmp)
cat("\n")


if(NVC>0){
cat("\n")
cat(" *** Variance-covariance matrix of the random-effects:\n" )
if(x$idiag==1){
Mat.cov <- diag(coef[(NEF+1):(NEF+NVC)])
colnames(Mat.cov) <-x$name.mat.cov 
rownames(Mat.cov) <-x$name.mat.cov 
Mat.cov[lower.tri(Mat.cov)] <- 0
Mat.cov[upper.tri(Mat.cov)] <- NA

print(Mat.cov,na.print="")
cat("\n")
}


if(x$idiag==0){
Mat.cov<-matrix(0,ncol=length(x$name.mat.cov),nrow=length(x$name.mat.cov))
colnames(Mat.cov) <-x$name.mat.cov 
rownames(Mat.cov) <-x$name.mat.cov 
Mat.cov[upper.tri(Mat.cov,diag=TRUE)]<-coef[(NEF+1):(NEF+NVC)]
Mat.cov <-t(Mat.cov)
Mat.cov[upper.tri(Mat.cov)] <- NA

print(Mat.cov,na.print="")
cat("\n")


}
}

std <- cbind(coef[NPM],se[NPM])
colnames(std) <-c("coef","se") 
rownames(std) <-" *** Residual standard error:"
if(NW>=1) {
nom <- paste(" *** Proportional variance coefficient",c(1:(x$ng-1)))
std <-cbind(coef[(NEF+NVC+1):NPM],se[(NEF+NVC+1):NPM]) 
rownames(std) <- c(nom,"Residual standard error")
colnames(std) <-c("coef","se") 
}
print(std)
cat("\n")

}

}

