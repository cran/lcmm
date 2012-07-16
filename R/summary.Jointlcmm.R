summary.Jointlcmm <-
function(object,...){
x <- object
if (!inherits(x, "Jointlcmm")) stop("use only with \"Jointlcmm\" objects")

#cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")



cat("Joint latent class model for quantitative longitudinal outcome", "\n")
cat("     fitted by maximum likelihood method", "\n")

cl <- x$call
cl$B <- NULL
cat(" \n")
dput(cl)
cat(" \n")

cat("Statistical Model:", "\n")
cat(paste("     Dataset:", x$dataset),"\n")
cat(paste("     Number of subjects:", x$specif[[2]]),"\n")
if(length(x$na.action))cat(paste("     Number of observations deleted:",length(x$na.action)),"\n")
cat(paste("     Number of repeated measures:", length(x$pred[,1])),"\n")
cat(paste("     Number of events: ", x$specif[[1]][7]," (",round(x$specif[[1]][7]*100/x$specif[[2]],2),"%)"),"\n")
cat(paste("     Number of latent classes:", x$specif[[3]]), "\n")
cat(paste("     Number of parameters:", length(x$best))," \n")


if (x$specif[[3]]>1) {
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
cat(     x$hazard[[3]]," \n")
}
if (x$hazard[[1]]==3)
{
cat("     M-splines constant baseline risk function with nodes \n")
cat(     x$hazard[[3]]," \n")
}



cat(" \n")
cat("Iteration process:", "\n")

if(x$conv==1) cat("     Convergence criteria satisfied")
if(x$conv==2) cat("     Maximum number of iteration reached without convergence")
if(x$conv==4|x$conv==12){

cat("     The program stopped abnormally. No results can be provided.\n")
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
cat(paste("     BIC:", round(-2*x$loglik+length(x$best)*log(x$specif[[2]]),2))," \n")
if(!is.na(x$CIstat)){
cat(paste("     Score test statistic for CI assumption: ", round(x$CIstat,3)," (p-value=",round((1-pchisq(x$CIstat,sum(x$specif[[4]]))),4),")" ,sep=""))
}
cat(" \n")
cat(" \n")
cat("Maximum Likelihood Estimates:", "\n")
cat(" \n")

NPROB <- x$specif[[1]][1]
nrisqtot <-x$specif[[1]][2]
nvarxevt <-x$specif[[1]][3]
NEF   <- x$specif[[1]][4]
NVC   <- x$specif[[1]][5]
NW    <- x$specif[[1]][6]
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

tmp <- cbind(round(coef[1:NPROB],5),round(se[1:NPROB],5),round(wald[1:NPROB],3),round(pwald[1:NPROB],5))
dimnames(tmp) <- list(names(coef)[1:NPROB], c("coef", "Se", "Wald", "p-value"))
cat("\n")
prmatrix(tmp,na.print="")
cat("\n")

}

if((nrisqtot+nvarxevt)>0){
cat(" *** Parameters in the proportional hazard model:\n" )

tmp <- cbind(round(coef[(NPROB+1):(NPROB+nrisqtot+nvarxevt)],5),round(se[(NPROB+1):(NPROB+nrisqtot+nvarxevt)],5),round(wald[(NPROB+1):(NPROB+nrisqtot+nvarxevt)],3),round(pwald[(NPROB+1):(NPROB+nrisqtot+nvarxevt)],5))
dimnames(tmp) <- list(names(coef)[(NPROB+1):(NPROB+nrisqtot+nvarxevt)], c("coef", "Se", "Wald", "p-value"))
cat("\n")
prmatrix(tmp,na.print="")
cat("\n")

}


cat(" *** Fixed effects in the longitudinal model:\n" )

tmp <- cbind(round(coef[(NPROB+nrisqtot+nvarxevt+1):(NEF)],5),round(se[(NPROB+nrisqtot+nvarxevt+1):(NEF)],5),round(wald[(NPROB+nrisqtot+nvarxevt+1):(NEF)],3),round(pwald[(NPROB+nrisqtot+nvarxevt+1):(NEF)],5))
dimnames(tmp) <- list(names(coef)[(NPROB+nrisqtot+nvarxevt+1):(NEF)], c("coef", "Se", "Wald", "p-value"))
cat("\n")
prmatrix(tmp,na.print="")
cat("\n")


if(NVC>0){
cat("\n")
cat(" *** Variance-covariance matrix of the random-effects:\n" )
if(x$specif[[8]]==1){
if (NVC>1) {
Mat.cov <- diag(coef[(NEF+1):(NEF+NVC)])
}else{
Mat.cov <- matrix(coef[(NEF+1)],ncol=1)
}
colnames(Mat.cov) <-x$Names[[2]][x$specif[[4]]!=0] 
rownames(Mat.cov) <-x$Names[[2]][x$specif[[4]]!=0] 
Mat.cov[lower.tri(Mat.cov)] <- 0
Mat.cov[upper.tri(Mat.cov)] <- NA

print(Mat.cov,na.print="")
cat("\n")
}


if(x$specif[[8]]==0){
Mat.cov<-matrix(0,ncol=sum(x$specif[[4]]),nrow=sum(x$specif[[4]]))
colnames(Mat.cov) <-x$Names[[2]][x$specif[[4]]!=0] 
rownames(Mat.cov) <-x$Names[[2]][x$specif[[4]]!=0] 
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
if((NW>=1)&(x$specif[[3]]>1)) {
nom <- paste("Proportional coefficient class",c(1:(x$specif[[3]]-1)),sep="")
std <-cbind(coef[(NEF+NVC+1):NPM],se[(NEF+NVC+1):NPM]) 
rownames(std) <- c(nom,"Residual standard error")
colnames(std) <-c("coef","se") 
}
print(std,na.print="")
cat("\n")

}

}

