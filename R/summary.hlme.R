summary.hlme <-
function(object,...){
x <- object
if (!inherits(x, "hlme")) stop("use only with \"hlme\" objects")

cat("Heterogenous linear mixed model", "\n")
cat("     fitted by maximum likelihood method", "\n")

cl <- x$call
cat(" \n")
dput(cl)
cat(" \n")

cat("Statistical Model:", "\n")
cat(paste("     Dataset:", x$dataset),"\n")
cat(paste("     Number of subjects:", x$ns),"\n")
cat(paste("     Number of observations:", length(x$pred[,1])),"\n")
cat(paste("     Number of latents classes:", x$ng), "\n")
cat(paste("     Number of parameters:", length(x$best))," \n")
cat(" \n")
cat("Iteration process:", "\n")

if(x$conv==1) cat("Convergence criteria satisfied")
if(x$conv==2) cat("Maximum number of iteration reached without convergence")
if(x$conv==4|x$conv==12) {
cat("The program stopped abnormally. No results can be displayed.\n")
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

cat("Maximum Likelihood Estimates:", "\n")
cat(" \n")

NPROB <- x$N[1]
NEF   <- x$N[2]
NVC   <- x$N[3]
NW    <- x$N[4]
NPM   <- length(x$best)


se <- rep(1,NPM)
if ((all.equal(x$conv,1)==T)==T){
#recuperation des indices de V
id <- 1:NPM
indice <- rep(id*(id+1)/2)
se <-sqrt(x$V[indice])
if (NVC>0) se[(NPROB+NEF+1):(NPROB+NEF+NVC)]<-1
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
cat("Fixed effects in the class-membership model:\n" )
cat("(the class of reference is the last class) \n")

tmp <- cbind(coef[1:NPROB],se[1:NPROB],wald[1:NPROB],pwald[1:NPROB])
dimnames(tmp) <- list(names(coef)[1:NPROB], c("coef", "Se", "Wald", "p-value"))
cat("\n")
prmatrix(tmp)
cat("\n")

}


cat("Fixed effects in the longitudinal model:\n" )

tmp <- cbind(coef[(NPROB+1):(NPROB+NEF)],se[(NPROB+1):(NPROB+NEF)],wald[(NPROB+1):(NPROB+NEF)],pwald[(NPROB+1):(NPROB+NEF)])
dimnames(tmp) <- list(names(coef)[(NPROB+1):(NPROB+NEF)], c("coef", "Se", "Wald", "p-value"))
cat("\n")
prmatrix(tmp)
cat("\n")


if(NVC>0){
cat("\n")
cat("Variance-covariance matrix of the random-effects:\n" )
if(x$idiag==1){
if (NVC>1) {
Mat.cov <- diag(coef[(NPROB+NEF+1):(NPROB+NEF+NVC)])
}else{
Mat.cov <- matrix(coef[(NPROB+NEF+1)],ncol=1)
}
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
Mat.cov[upper.tri(Mat.cov,diag=TRUE)]<-coef[(NPROB+NEF+1):(NPROB+NEF+NVC)]
Mat.cov <-t(Mat.cov)
Mat.cov[upper.tri(Mat.cov)] <- NA

print(Mat.cov,na.print="")
cat("\n")
}
}

std <- cbind(coef[NPM],se[NPM])
colnames(std) <-c("coef","se") 
rownames(std) <-"Residual standard error:"
if(NW>=1) {
nom <- paste("Proportional variance coefficient",c(1:(x$ng-1)))
std <-cbind(coef[(NPROB+NEF+NVC+1):NPM],se[(NPROB+NEF+NVC+1):NPM]) 
rownames(std) <- c(nom,"Residual standard error")
colnames(std) <-c("coef","se") 
}
print(std)
cat("\n")
}
}


