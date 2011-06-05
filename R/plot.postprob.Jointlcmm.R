plot.postprob.Jointlcmm <-
function(x,...){
if (!inherits(x, "Jointlcmm")) stop("use only with \"Jointlcmm\" objects")
if(x$conv==1|x$conv==2) {
   if(sum(is.na(x$pprob[,3:(x$ng+2)]))==0){
      if(x$ng>1){
        for(i in 1 : x$ng){
          dev.new()
          xlab1 <- paste("class",i)
          title.hist <- paste("Posterior probabilities given longitudinal and time-to-event data in class",i)
          hist(x$pprob[,i+2],prob=TRUE,xlab=xlab1,main=title.hist,...)
        }
      }else{stop("plot only for ng > 1")}

   }else{
      stop("Error in the computation of posterior class-membership probabilities given all the information")
   }       
   }else{
      cat("Output can not be produced since the program stopped abnormally.")
      stop("Pease check the data & model specification")
   }
}



plot.postprob <- function(x,...) UseMethod("plot.postprob")
