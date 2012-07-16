plot.lcmm <-
function(x,...){
if (!inherits(x, "lcmm")) stop("use only with \"lcmm\" objects")


if (x$linktype==3) 
{
cat("Residual and prediction plots are not available yet for threshold mixed models \n")
}else{
if(x$conv==1|x$conv==2) {
dev.new()
plot(x$pred[,3]~x$pred[,2],xlab="marginal predictions",ylab="marginal residuals",main="marginal residuals versus marginal predictions",...)
dev.new()
plot(x$pred[,5]~x$pred[,4],xlab="subject-specific predictions",ylab="subject-specific residuals",main="subject-specific residuals versus subject-specific predictions",...)
}else{
cat("Output can not be produced since the program stopped abnormally.")
}
}
}
