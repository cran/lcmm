plot.hlme <-
function(x,...){
if (!inherits(x, "hlme")) stop("use only with \"hlme\" objects")
if(x$conv==1|x$conv==2) {
dev.new()
plot(x$pred[,3]~x$pred[,2],xlab="marginal predictions",ylab="marginal residuals",main="marginal residuals versus marginal predictions",...)
dev.new()
plot(x$pred[,5]~x$pred[,4],xlab="subject-specific predictions",ylab="subject-specific residuals",main="subject-specific residuals versus subject-specific predictions",...)
}else{
cat("Output can not be produced since the program stopped abnormally.")
stop("Pease check the data & model specification")
}
}

