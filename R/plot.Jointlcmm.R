plot.Jointlcmm <- function(x,...)
{
 if (!inherits(x, "Jointlcmm")) stop("use only with \"Jointlcmm\" objects")
 if(x$conv==1|x$conv==2)
 {
  dots <- list(...)
  
   if(length(list(...)$main)) 
   {
    title1 <- as.character(eval(match.call()$main))
    dots <- dots[setdiff(names(dots),"main")]
   }
   else title1 <- c("marginal residuals versus marginal predictions", "subject-specific residuals versus subject-specific predictions")
   title1 <- rep(title1,length.out=2)                         
   
   if(length(list(...)$xlab)) 
   {
    xlab1 <- as.character(eval(match.call()$xlab))
    dots <- dots[setdiff(names(dots),"xlab")]
   }
   else xlab1 <- c("marginal predictions","subject-specific predictions")
   xlab1 <- rep(xlab1,length.out=2)

   if(length(list(...)$ylab)) 
   {
    ylab1 <- as.character(eval(match.call()$ylab))
    dots <- dots[setdiff(names(dots),"ylab")]
   }
   else ylab1 <- c("marginal residuals","subject-specific residuals")
   ylab1 <- rep(ylab1,length.out=2)    
   

  names.plot <- c("adj","ann","asp","axes","bg","bty","cex","cex.axis","cex.lab","cex.main","cex.sub","col","col.axis",
  "col.lab","col.main","col.sub","crt","err","family","fig","fin","font","font.axis","font.lab","font.main","font.sub",
  "frame.plot","lab","las","lend","lheight","ljoin","lmitre","lty","lwd","mai","main","mar","mex","mgp","mkh","oma",
  "omd","omi","pch","pin","plt","ps","pty","smo","srt","sub","tck","tcl","type","usr","xaxp","xaxs","xaxt","xlab",
  "xlim","xpd","yaxp","yaxs","yaxt","ylab","ylbias","ylim") 
  dots.plot <- dots[intersect(names(dots),names.plot)]  
  
  par(mfrow=c(1,2))
  do.call("plot",c(dots.plot,list(x=x$pred[,2],y=x$pred[,3],main=title1[1],xlab=xlab1[1],ylab=ylab1[1])))
  do.call("plot",c(dots.plot,list(x=x$pred[,4],y=x$pred[,5],main=title1[2],xlab=xlab1[2],ylab=ylab1[2])))
 }
 else
 {
  cat("Output can not be produced since the program stopped abnormally.")
  stop("Pease check the data & model specification")
 }
}

