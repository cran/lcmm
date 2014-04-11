

plot.baselinerisk.Jointlcmm <- function(x,legend.loc="topleft",legend,add=FALSE,...)
{
  if(missing(x)) stop("The argument x should be specified")
  if (!inherits(x, "Jointlcmm")) stop("use only with \"Jointlcmm\" objects")
  if(is.na(as.logical(add))) stop("add should be TRUE or FALSE")
  
  if((x$conv==1|x$conv==2)& (sum(is.na(x$predSurv)==0)))
  {
   ng <- x$specif[[3]]
  
   dots <- list(...)

   if(length(list(...)$main)) 
   {
    title1 <- as.character(eval(match.call()$main))
    dots <- dots[setdiff(names(dots),"main")]
   }
   else title1 <- "Class-specific baseline risk functions"

   if(length(list(...)$col)) 
   {
    color <- as.vector(eval(match.call()$col))
    dots <- dots[-which(names(dots)=="col")]
   }
   else  color <- 1:ng                            
 
   if(length(list(...)$type))    
   {
    type1 <- eval(match.call()$type)
    dots <- dots[-which(names(dots)=="type")]
   }
   else  type1 <- "l"
      
   if(length(list(...)$xlab)) 
   {
    xlab1 <- as.character(eval(match.call()$xlab))
    dots <- dots[setdiff(names(dots),"xlab")]
   }
   else xlab1 <- "Time"
   
   if(length(list(...)$ylab)) 
   {
    ylab1 <- as.character(eval(match.call()$ylab))
    dots <- dots[setdiff(names(dots),"ylab")]
   }
   else ylab1 <- "Baseline risk function" 
   
   if(length(list(...)$lty))    
   {
    lty1 <- eval(match.call()$lty)
    dots <- dots[-which(names(dots)=="lty")]
   }
   else  lty1 <- 1:ng  
   
#   if("legend" %in% names(match.call())) 
#   {   
#    nomsleg <- eval(match.call()$legend)
#    dots <- dots[setdiff(names(dots),"legend")]
#   }
#   else nomsleg <- paste("class",1:ng,sep="")
   if(missing(legend)) legend <- paste("class",1:ng,sep="")
   
   if(length(list(...)$box.lty)) 
   {
    box.lty1 <- as.integer(eval(match.call()$box.lty))
    dots <- dots[setdiff(names(dots),"box.lty")]
   }
   else box.lty1 <- 0
   
   if(length(list(...)$inset)) 
   {
    inset1 <- eval(match.call()$inset)
    dots <- dots[setdiff(names(dots),"inset")]
   }
   else inset1 <- c(0.05,0.05)
   

   
  names.plot <- c("adj","ann","asp","axes","bg","bty","cex","cex.axis","cex.lab","cex.main","cex.sub","col","col.axis",
  "col.lab","col.main","col.sub","crt","err","family","fig","fin","font","font.axis","font.lab","font.main","font.sub",
  "frame.plot","lab","las","lend","lheight","ljoin","lmitre","lty","lwd","mai","main","mar","mex","mgp","mkh","oma",
  "omd","omi","pch","pin","plt","ps","pty","smo","srt","sub","tck","tcl","type","usr","xaxp","xaxs","xaxt","xlab",
  "xlim","xpd","yaxp","yaxs","yaxt","ylab","ylbias","ylim") 
  dots.plot <- dots[intersect(names(dots),names.plot)]
  
  if(!isTRUE(add))
  {
   do.call("matplot",c(dots.plot,list(x=x$predSurv[,1],y=x$predSurv[,(1+1:x$specif[[3]])],xlab=xlab1,ylab=ylab1,main=title1,type=type1,col=color)))
  }
  else
  {
   do.call("matlines",c(dots.plot,list(x=x$predSurv[,1],y=x$predSurv[,(1+1:x$specif[[3]])],type=type1,col=color)))
  }
  
  names.legend <- c("fill","border","lty","lwd","pch","angle","density","bg","box.lwd",   
  "box.lty","box.col","pt.bg","cex","pt.cex","pt.lwd","xjust","yjust","x.intersp","y.intersp","adj","text.width",
  "text.col","text.font","merge","trace","plot","ncol","horiz","title","xpd","title.col","title.adj","seg.len")    
  dots.leg <- dots[intersect(names(dots),names.legend)]
  
#   for (g in 1:x$specif[[3]])
#   {
#    Y <- x$predSurv[,(1+g)]
#    if(g==1) plot(Y~x$predSurv[,1],col=color[g],type=type1,xlim=xlim1,ylim=ylim1,xlab=xlab1,ylab=ylab1,pch=pch1[g],bg=bg1[g],main=title1,lty=lty1[g],lwd=lwd1[g],frame.plot=frame.plot1,mgp=mgp1,axes=axes1,xaxt=xaxt1,yaxt=yaxt1)
#    else lines(Y~x$predSurv[,1],col=color[g],type=type1,pch=pch1[g],lty=lty1[g],lwd=lwd1[g],bg=bg1[g])
#   }
                                    
  	if(!is.null(legend))
    {
     #if(all(type1=="p") & !all(is.na(pch1))) lty1 <- 0
     do.call("legend",c(dots.leg,list(x=legend.loc,legend=legend,col=color,box.lty=box.lty1,inset=inset1,lty=lty1)))
    }
  }
  else
  {
   cat("Output can not be produced. The program stopped abnormally or there was an error in the computation of the estimated baseline risk functions and survival functions.\n")
  }
}


plot.baselinerisk <- function(x,legend.loc="topleft",legend,add=FALSE,...) UseMethod("plot.baselinerisk")

