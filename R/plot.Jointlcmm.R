plot.Jointlcmm <-
function(x,...){
if (!inherits(x, "Jointlcmm")) stop("use only with \"Jointlcmm\" objects")
if(x$conv==1|x$conv==2) {
dev.new()
plot(x$pred[,2]~x$pred[,1],xlab="marginal predictions",ylab="marginal residuals",main="marginal residuals versus marginal predictions",...)
dev.new()
plot(x$pred[,4]~x$pred[,3],xlab="subject-specific predictions",ylab="subject-specific residuals",main="subject-specific residuals versus subject-specific predictions",...)
}else{
cat("Output can not be produced since the program stopped abnormally.")
stop("Pease check the data & model specification")
}
}

