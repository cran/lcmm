plot.hlme <-
function(x,...){
if (!inherits(x, "hlme")) stop("use only with \"hlme\" objects")
dev.new()
plot(x$pred[,2]~x$pred[,1],xlab="marginal predictions",ylab="marginal residuals",main="marginal residuals versus marginal predictions",...)
dev.new()
plot(x$pred[,4]~x$pred[,3],xlab="subject-specific predictions",ylab="subject-specific residuals",main="subject-specific residuals versus subject-specific predictions",...)
}

