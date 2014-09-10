plot.Diffepoce <- function(x,...)
{
    if (!inherits(x, "Diffepoce")) stop("use only with \"Diffepoce\" objects")

   if (all(is.na(x$DiffEPOCE[,3]))|all(is.na(x$DiffEPOCE[,4]))) stop("can't produce the plot with missing differences in EPOCE")
    
if (all(is.infinite(x$DiffEPOCE[,3]))|all(is.infinite(x$DiffEPOCE[,4]))) stop("can't produce the plot with infinite differences in EPOCE")
    
   dots <- list(...)

   if(length(list(...)$main)) 
   {
    title1 <- as.character(eval(match.call()$main))
    dots <- dots[setdiff(names(dots),"main")]
   }
   else title1 <- "Difference in EPOCE estimates"  
   
   if(length(list(...)$col)) 
   {
    color <- as.vector(eval(match.call()$col))
    dots <- dots[-which(names(dots)=="col")]
   }
   else  color <- c("black","black","black","lightgrey")                             

   if(length(list(...)$type))    
   {
    type1 <- eval(match.call()$type)
    dots <- dots[-which(names(dots)=="type")]
   }
   else  type1 <- c("b","b","b","l")

   if(length(list(...)$pch))    
   {
    pch1 <- eval(match.call()$pch)
    dots <- dots[-which(names(dots)=="pch")]
   }
   else  pch1 <- c(18,18,18,NA)

   if(length(list(...)$lty))    
   {
    lty1 <- eval(match.call()$lty)
    dots <- dots[-which(names(dots)=="lty")]
   }
   else  lty1 <- c(1,3,3,1)

    if(length(list(...)$ylim)) 
        {
            ylim1 <- eval(match.call()$ylim)
            dots <- dots[setdiff(names(dots),"ylim")]
        }
    else
        {
            if(all(is.na(x$DiffEPOCE[,2:4])) | all(is.infinite(x$DiffEPOCE[,2:4])))
                {
                    ylim1 <- c(-1000,1000)
                }
            else
                {
                    ylim1 <- c(min(x$DiffEPOCE[,2:4][!(is.na(x$DiffEPOCE[,2:4])) & is.finite(x$DiffEPOCE[,2:4])]),max(x$DiffEPOCE[,2:4][!(is.na(x$DiffEPOCE[,2:4])) & is.finite(x$DiffEPOCE[,2:4])]))
                }
        }
   
   if(length(list(...)$xlab)) 
   {
    xlab1 <- as.character(eval(match.call()$xlab))
    dots <- dots[setdiff(names(dots),"xlab")]
   }
   else xlab1 <- "prediction time"

   if(length(list(...)$ylab)) 
   {
    ylab1 <- as.character(eval(match.call()$ylab))
    dots <- dots[setdiff(names(dots),"ylab")]
   }
   else 
   {
  	if(x$new.data==FALSE) 
    {
		 ylab1 <- expression(Delta(CVPOL))
	  }
    else
    {
		 ylab1 <- expression(Delta(MPOL))	
    }
   }   


   
  names.plot <- c("adj","ann","asp","axes","bg","bty","cex","cex.axis","cex.lab","cex.main","cex.sub","col","col.axis",
  "col.lab","col.main","col.sub","crt","err","family","fig","fin","font","font.axis","font.lab","font.main","font.sub",
  "frame.plot","lab","las","lend","lheight","ljoin","lmitre","lty","lwd","mai","main","mar","mex","mgp","mkh","oma",
  "omd","omi","pch","pin","plt","ps","pty","smo","srt","sub","tck","tcl","type","usr","xaxp","xaxs","xaxt","xlab",
  "xlim","xpd","yaxp","yaxs","yaxt","ylab","ylbias","ylim") 
  dots.plot <- dots[intersect(names(dots),names.plot)]
  
  do.call("matplot",c(dots.plot,list(x=x$DiffEPOCE[,1],y=cbind(x$DiffEPOCE[,2:4],0),xlab=xlab1,ylab=ylab1,main=title1,type=type1,col=color,pch=pch1,lty=lty1,ylim=ylim1)))
}



