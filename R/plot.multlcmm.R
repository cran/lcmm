plot.multlcmm <- function(x,...)
{
  if (!inherits(x, "multlcmm")) stop("use only with \"multlcmm\" objects")

  if (any(x$linktype==3))
  {
   cat("Residual and prediction plots are not available yet for threshold mixed models \n")
  }
  else
  {
   if(x$conv==1|x$conv==2)
   {
    dev.new()
    plot(x$pred[,"resid_m"]~x$pred[,"pred_m"],xlab="marginal predictions",ylab="marginal residuals",main="marginal residuals versus marginal predictions",...)
    dev.new()
    plot(x$pred[,"resid_ss"]~x$pred[,"pred_ss"],xlab="subject-specific predictions",ylab="subject-specific residuals",main="subject-specific residuals versus subject-specific predictions",...)
   }
   else
   {
    cat("Output can not be produced since the program stopped abnormally.")
   }
  }
}