plot.dynpred <- function(x,subject=NULL,landmark=NULL,horizon=NULL,add=FALSE,...)
{
 if(missing(x)) stop("The argument \'x\' is missing.")
 if(!inherits(x,"dynpred")) stop("use only with \'dynpred\' object")
 if(length(x)==1){if(is.na(x)) stop("x is NA")}
  
 
 if(!is.null(landmark))
 {
  if(!all(landmark %in% x$pred[,2])) stop(paste("Dynamic predictions have only been calculated on landmark ",paste(unique(x$pred[,2]),collapse=","),sep=""))
 }
 
 if(!is.null(horizon))
 {
  if(!all(horizon %in% x$pred[,3])) stop(paste("Dynamic predictions have only been calculated on horizon ",paste(unique(x$pred[,3]),collapse=","),sep=""))
 }
 
 if(!is.null(subject))
 {
  if(!all(subject %in% x$pred[,1])) stop(paste("Dynamic predictions have only been calculated for subject ",paste(unique(x$pred[,1]),collapse=","),sep="")) 
 }
 
 if(length(landmark)>1 & length(horizon)>1) stop("Either landmark or horizon should be of length 1")
 
 if(is.null(landmark) & is.null(horizon))
 {
  if(length(unique(x$pred[,2]))>1 & length(unique(x$pred[,3]))>1) stop("Only predictions for a fixed landmark or horizon can be plotted. Please specify it in \'landmark\' or \'horizon\'.")
 }
 
 if(is.null(landmark)) landmark <- unique(x$pred[,2])
 if(is.null(horizon)) horizon <- unique(x$pred[,3])
 if(is.null(subject)) subject <- unique(x$pred[,1])
 
 
 
 ##parametres graphiques
 dots <- list(...)

 names.plot <- c("adj","ann","asp","axes","bg","bty","cex","cex.axis","cex.lab","cex.main","cex.sub","col.axis",
 "col.lab","col.main","col.sub","crt","err","family","fig","fin","font","font.axis","font.lab","font.main","font.sub",
 "frame.plot","lab","las","lend","lheight","ljoin","lmitre","lwd","mai","main","mar","mex","mgp","mkh","oma",
 "omd","omi","pin","plt","ps","pty","smo","srt","sub","tck","tcl","usr","xaxp","xaxs","xaxt","xlab",
 "xlim","xpd","yaxp","yaxs","yaxt","ylab","ylbias","ylim") #col,pch,type,lty

 names.legend <- c("adj","angle","bg","border","box.col","box.lty","box.lwd", 
 "cex","col","density","fill","horiz","lty","lwd","merge", 
 "ncol","pch","plot","pt.bg","pt.cex","pt.lwd","seg.len", 
 "text.col","text.font","text.width","title","title.adj", 
 "title.col","trace","x.intersp","xjust","xpd","y.intersp","yjust")


 if(!length(dots$xlab))
 {
  dots$xlab <- "time"
 }

 if(!length(dots$ylab))
 {
  dots$ylab <- "Longitudinal marker"
 }   

  
 if(length(dots$col))
 {
  if(length(dots$col) %in% c(length(subject),length(subject)*2))
  {
   col1 <- rep(dots$col,length.out=2*length(subject))
  }
  else
  {
   col1 <- rep(rainbow(length(subject)),2)
  }
 }
 else
 {
  col1 <- rep(rainbow(length(subject)),2) 
 }   
 
 if(length(dots$pch))
 {
  if(length(dots$pch) %in% c(length(subject),length(subject)*2))
  {
   pch1 <- rep(dots$pch,length.out=2*length(subject))
  }
  else
  {
   pch1 <- c(4+(0:(length(subject)-1)*2),3+(0:(length(subject)-1)*2)) 
  }
 }
 else
 {
  pch1 <- c(4+(0:(length(subject)-1)*2),3+(0:(length(subject)-1)*2)) 
 }
 
 if(!length(dots$lwd))
 {
  dots$lwd <- 1
 }
  
 dots.plot <- dots[intersect(names(dots),names.plot)]  
 dots.leg <- dots[intersect(names(dots),names.legend)]

 if(length(dots$ylab)==2)
 {
  dots.plot$ylab <- dots$ylab[1] 
 }   
 
 if(ncol(x$pred)==6)
 {
  bout <- min(abs(diff(x$newdata[,"time"])))/10*ifelse(all(dots.plot$lwd<3),max(dots.plot$lwd),2)
 }

 oldmar <- par()$mar
 par(mar=c(5.1,4.1,4.1,4.1))
 #on.exit(par(mar=oldmar))  
  
 if(length(horizon)==1)
 {
  res <- x$pred[which((x$pred[,2] %in% landmark) & (x$pred[,3]==horizon)),,drop=FALSE]
  newdata <- x$newdata[which(x$newdata[,"time"]<max(landmark)),,drop=FALSE]
   
  if(all(is.na(res[,4]))) stop("There is no prediction in 'x'")
 
  if(!length(dots.plot$xlim))
  {
   dots.plot$xlim <- range(c(newdata$time[which(newdata$id %in% subject)],landmark))
  }

  if(!length(dots.plot$ylim))
  {
   dots.plot$ylim <- range(newdata$y[which(newdata$id %in% subject)])
  } 
  if(dots.plot$ylim[1]==dots.plot$ylim[2])
  {
   dots.plot$ylim <- c(dots.plot$ylim[1]-dots.plot$ylim[1]/2,dots.plot$ylim[1]+dots.plot$ylim[1]/2) 
  }    
   
  if(!length(dots.plot$main))
  {
   dots.plot$main <- paste("Predictions for horizon",horizon) 
  }  
 
  fromptoy <- function(p)
  {
   return((dots.plot$ylim[2]-dots.plot$ylim[1])*p+dots.plot$ylim[1]) 
  }    
  
  iok <- 0
   
  for(i in 1:length(subject))
  {
   #selectionner les donnees de ce sujet
   newdatai <- newdata[which(newdata$id==subject[i]),,drop=FALSE]
   predi <- res[which(res[,1]==subject[i]),,drop=FALSE] 
   
   if(all(is.na(predi[,4])))
   {
    next
   }
   else
   {
    iok <- iok+1 
   }
    
   #tracer les observations
   if((i==1 | iok==1) & !isTRUE(add))
   {
    do.call(plot,c(dots.plot,list(x=newdatai$time,y=newdatai$y,type="p",col=col1[i],pch=pch1[i])))
   }
   else
   {
    do.call(points,c(dots.plot[setdiff(names(dots.plot),"axes")],list(x=newdatai$time,y=newdatai$y,col=col1[i],pch=pch1[i])))
   }
   
   #tracer les predictions
    do.call(points,c(dots.plot[setdiff(names(dots.plot),"axes")],list(x=predi[,2],y=fromptoy(predi[,4]),col=col1[length(subject)+i],pch=pch1[length(subject)+i])))
    
   #tracer les IC
   if(ncol(x$pred)==6)
   {
    do.call(segments,c(dots.plot[setdiff(names(dots.plot),c("type","axes"))],
    list(x0=predi[,2,drop=FALSE],
         y0=fromptoy(predi[,5,drop=FALSE]),
         x=predi[,2,drop=FALSE],
         y=fromptoy(predi[,6,drop=FALSE]),
         col=col1[length(subject)+i])))
   #tracer les bouts des IC :
   do.call(segments,c(dots.plot[setdiff(names(dots.plot),c("type","axes"))],
   list(x0=predi[,2,drop=FALSE]-bout,
        y0=fromptoy(predi[,5,drop=FALSE]),
        x=predi[,2,drop=FALSE]+bout,
        y=fromptoy(predi[,5,drop=FALSE]),
        col=col1[length(subject)+i])))
   do.call(segments,c(dots.plot[setdiff(names(dots.plot),c("type","axes"))],
   list(x0=predi[,2,drop=FALSE]-bout,
        y0=fromptoy(predi[,6,drop=FALSE]),
        x=predi[,2,drop=FALSE]+bout,
        y=fromptoy(predi[,6,drop=FALSE]),
        col=col1[length(subject)+i])))
   }    
  }
   
  #tracer l'axe pour les probas   
  if(isTRUE(dots$axes) | (length(dots$axes)==0)) 
  {
   if((length(dots$yaxt)==0) | (length(dots$yaxt) & isTRUE(dots$yaxt[length(dots$yaxt)]!="n")))
   {   
    axis(side=4,at=fromptoy(pretty(seq(0,1))),labels=pretty(seq(0,1)))  
   }
  } 
  if(length(dots$ylab)==2)
  {
   mtext(side=4,text=dots$ylab[2],line=2) 
  }     
  else
  {  
   mtext(side=4,text="Probability of event",line=2)  
  }  
 }
 else #cas d'un seul landmark et plusieurs horizons
 {
  res <- x$pred[which((x$pred[,2]==landmark) & (x$pred[,3] %in% horizon)),,drop=FALSE]
  newdata <- x$newdata[which(x$newdata[,"time"]<max(landmark)),,drop=FALSE] 
   
  if(all(is.na(res[,4]))) stop("There is no prediction in 'x'") 

  if(!length(dots.plot$xlim))
  {
   dots.plot$xlim <- range(c(newdata$time[which(newdata$id %in% subject)],landmark+max(res[,3])))
  }
   
  if(!length(dots.plot$ylim))
  {
   dots.plot$ylim <- range(newdata$y[which(newdata$id %in% subject)])
  } 
  if(dots.plot$ylim[1]==dots.plot$ylim[2])
  {
   dots.plot$ylim <- c(dots.plot$ylim[1]-dots.plot$ylim[1]/2,dots.plot$ylim[1]+dots.plot$ylim[1]/2) 
  }  
   
  if(!length(dots.plot$main))
  {
   dots.plot$main <- paste("Predictions at landmark",landmark) 
  }  
   
  fromptoy <- function(p)
  {
   return((dots.plot$ylim[2]-dots.plot$ylim[1])*p+dots.plot$ylim[1]) 
  }  

  iok <- 0
   
  for(i in 1:length(subject))
  {
   newdatai <- newdata[which(newdata$id==subject[i]),,drop=FALSE]
   predi <- res[which(res[,1]==subject[i]),,drop=FALSE]
   
   if(all(is.na(predi[,4])))
   {
    next
   }
   else
   {
    iok <- iok+1 
   }
   
   #tracer les observations
   if((i==1 | iok==1) & !isTRUE(add))
   {
    do.call(plot,c(dots.plot,list(x=newdatai$time,y=newdatai$y,type="p",col=col1[i],pch=pch1[i])))
   }
   else
   {
    do.call(points,c(dots.plot[setdiff(names(dots.plot),"axes")],list(x=newdatai$time,y=newdatai$y,col=col1[i],pch=pch1[i])))
   }
   
   #tracer les predictions
   do.call(lines,c(dots.plot[setdiff(names(dots.plot),"axes")],list(x=predi[,2]+predi[,3],y=fromptoy(predi[,4]),col=col1[length(subject)+i])))

   #tracer les IC
   if(ncol(x$pred)==6)
   {
    do.call(matlines,c(dots.plot[setdiff(names(dots.plot),"axes")],
    list(x=predi[,2]+predi[,3],y=cbind(fromptoy(predi[,5]),fromptoy(predi[,6])),col=col1[length(subject)+i],lty=2)))
   }
  }
   
  #tracer l'axe pour les probas 
  if(isTRUE(dots$axes) | (length(dots$axes)==0)) 
  {
   if((length(dots$yaxt)==0) | (length(dots$yaxt) & isTRUE(dots$yaxt[length(dots$yaxt)]!="n")))
   {   
    axis(side=4,at=fromptoy(pretty(seq(0,1))),labels=pretty(seq(0,1)))  
   }
  }  
  if(length(dots$ylab)==2)
  {
   mtext(side=4,text=dots$ylab[2],line=2) 
  }     
  else
  {  
   mtext(side=4,text="Probability of event",line=2)  
  } 
 }
}

