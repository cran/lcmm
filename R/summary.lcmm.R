summary.lcmm <-
function(object,...){
x <- object
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
ntrtot <- 2
cat("     Link function: linear"," \n")
}
if (x$linktype==1)
{
ntrtot <- 4
cat("     Link function: Standardised Beta CdF"," \n")
}
if (x$linktype==2) {
ntrtot <- length(x$linknodes)+2
cat("     Link function: Quadratic I-splines with nodes"," \n")
#cat(paste("      ",x$linknodes)," \n")
cat(     x$linknodes," \n")
}


if (x$linktype==3) {
ntrtot <- sum(x$ide==1)
cat("     Link function: thresholds"," \n")
}

cat(" \n")
cat("Iteration process:", "\n")

if(x$conv==1) cat("     Convergence criteria satisfied")
if(x$conv==2) cat("     Maximum number of iteration reached without convergence")
if(x$conv==4|x$conv==12) {
cat("     The program stopped abnormally. No results can be displayed.\n")
}
else{

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
ncor <- 0
if(length(x$N)>5)
{
 ncor <- x$N[6]
}
if (x$Ydiscrete==1&ncor==0){
cat(paste("     Discrete posterior log-likelihood:", round(x$discrete_loglik,2))," \n")
cat(paste("     Discrete AIC:", round(-2*x$discrete_loglik+2*length(x$best),2))," \n")
cat(" \n")
cat(paste("     Mean discrete AIC per subject:",round((-x$discrete_loglik+length(x$best))/as.double(x$ns),4))," \n")
cat(paste("     Mean UACV per subject:",round(x$UACV,4))," \n")
cat(paste("     Mean discrete LL per subject:",round(x$discrete_loglik/as.double(x$ns),4))," \n")
}
cat(" \n")

cat("Maximum Likelihood Estimates:", "\n")
cat(" \n")

NPROB <- x$N[1]
NEF   <- x$N[2]
NVC   <- x$N[3]
NW    <- x$N[4]
NPM   <- length(x$best)


se <- rep(NA,NPM)
if ((all.equal(x$conv,1)==T)==T){
#recuperation des indices de V
id <- 1:NPM
indice <- rep(id*(id+1)/2)
se <-sqrt(x$V[indice])
if (NVC>0) se[(NPROB+NEF+1):(NPROB+NEF+NVC)]<-NA
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

tmp <- cbind(round(coef[1:NPROB],5),round(se[1:NPROB],5),round(wald[1:NPROB],3),round(pwald[1:NPROB],5))
dimnames(tmp) <- list(names(coef)[1:NPROB], c("coef", "Se", "Wald", "p-value"))
cat("\n")
prmatrix(tmp,na.print="")
cat("\n")
}



cat("Fixed effects in the longitudinal model:\n" )

tmp <- cbind(round(coef[(NPROB+1):(NPROB+NEF)],5),round(se[(NPROB+1):(NPROB+NEF)],5),round(wald[(NPROB+1):(NPROB+NEF)],3),round(pwald[(NPROB+1):(NPROB+NEF)],5))
tmp <- rbind(c(0,NA,NA,NA),tmp)
interc <- "intercept"
if (x$ng>1){
interc <- paste(interc,"class1")
}
interc <- paste(interc,"(not estimated)")
dimnames(tmp) <- list(c(interc,names(coef)[(NPROB+1):(NPROB+NEF)]), c("coef", "Se", "Wald", "p-value"))
cat("\n")

tTable <- tmp

prmatrix(tmp,na.print="")
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
colnames(Mat.cov) <-x$Xnames[x$idea0==1]
rownames(Mat.cov) <-x$Xnames[x$idea0==1]
Mat.cov[lower.tri(Mat.cov)] <- 0
Mat.cov[upper.tri(Mat.cov)] <- NA
print(Mat.cov,na.print="")
cat("\n")
}
if(x$idiag==0){
Mat.cov<-matrix(0,ncol=sum(x$idea0),nrow=sum(x$idea0))
colnames(Mat.cov) <-x$Xnames[x$idea0==1]
rownames(Mat.cov) <-x$Xnames[x$idea0==1]
Mat.cov[upper.tri(Mat.cov,diag=TRUE)]<-coef[(NPROB+NEF+1):(NPROB+NEF+NVC)]
Mat.cov <-t(Mat.cov)
Mat.cov[upper.tri(Mat.cov)] <- NA
print(Mat.cov,na.print="")
cat("\n")
}
}

std <- NULL
nom <- NULL
if(NW>=1) {
nom <- paste("Proportional coefficient class",c(1:(x$ng-1)),sep="")
std <-cbind(coef[(NPROB+NEF+NVC+1):(NPROB+NEF+NVC+NW)],se[(NPROB+NEF+NVC+1):(NPROB+NEF+NVC+NW)])
}
if(ncor==2) {
nom <- c(nom,"AR correlation parameter:","AR standard error:")
std <-rbind(std,c(coef[(NPROB+NEF+NVC+NW+ntrtot+1)],se[(NPROB+NEF+NVC+NW+ntrtot+1)]),c(abs(coef[(NPROB+NEF+NVC+NW+ntrtot+2)]),se[(NPROB+NEF+NVC+NW+ntrtot+2)]))
}
if(ncor==1) {
nom <- c(nom,"BM standard error:")
std <-rbind(std,c(abs(coef[(NPROB+NEF+NVC+NW+ntrtot+1)]),se[(NPROB+NEF+NVC+NW+ntrtot+1)]))
}
if (!is.null(std)) { 
rownames(std) <- nom
colnames(std) <-c("coef","se") 
print(std, na.print="")
cat("\n")
}

cat("Residual standard error (not estimated) = 1\n")
cat("\n")

cat("Parameters of the link function:\n" )
if (x$linktype==3 & ntrtot != (x$linknodes[2]-x$linknodes[1])) {
temp <- (x$linknodes[1]:(x$linknodes[2]-1))*(1-x$ide)
cat("(the following levels are not observed in the data: ",temp[temp!=0],"\n")
cat("so that the number of parameters in the threshold transformation is reduced to",ntrtot,") \n")
}

tmp <- cbind(round(coef[(NPM-ntrtot+1-ncor):(NPM-ncor)],5),round(se[(NPM-ntrtot+1-ncor):(NPM-ncor)],5),round(wald[(NPM-ntrtot+1-ncor):(NPM-ncor)],3),round(pwald[(NPM-ntrtot+1-ncor):(NPM-ncor)],5))
dimnames(tmp) <- list(names(coef)[(NPM-ntrtot+1-ncor):(NPM-ncor)], c("coef", "Se", "Wald", "p-value"))
cat("\n")
prmatrix(tmp,na.print="")
cat("\n")

return(invisible(tTable))
}
}
