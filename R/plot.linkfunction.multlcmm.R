
plot.linkfunction.multlcmm <- function(x,legend.loc="topleft",...)
{
  if(missing(x)) stop("The argument x should be specified")
  if (!inherits(x, "multlcmm")) stop("use only with \"multlcmm\" objects")
  if(x$conv %in% c(1,2))
  {
   ny <- length(x$Ynames)
   limx <- c(min(x$estimlink[,2*1:ny]),max(x$estimlink[,2*1:ny]))
   nsim <- length(x$estimlink[,1])
   
   if(length(list(...)$main)) {title1 <- as.character(eval(match.call()$title))}
   else title1 <- "Estimated link functions"
   
   if(length(list(...)$col)) {color <- as.vector(eval(match.call()$col))}
   else  color <- rainbow(ny)
     
   loc.grad <- function(y.grad,yk)
   {
    (nsim*(y.grad-min(x$linknodes[,yk]))-y.grad+max(x$linknodes[,yk]))/(max(x$linknodes[,yk])-min(x$linknodes[,yk]))
   }
   
   dev.new()
   par(mar=c(5,ny+1,4,ny+1)+0.2)
   
   plot(x=x$estimlink[,2],y=1:nsim,type="l",xlim=limx,ylim=c(1,nsim),xlab="Latent process",ylab="",main=title1,yaxt="n",col=color[1],frame.plot=FALSE)
   axis(1)
   y.grad <- pretty(min(x$linknodes[,1]):max(x$linknodes[,1]))
   y.grad[1] <- round(min(x$linknodes[,1]),2)
   y.grad[length(y.grad)] <- round(max(x$linknodes[,1]),2)
   axis(2,at=loc.grad(y.grad,1),labels=y.grad,col=color[1], col.axis=color[1], cex.axis=0.8)
   for (i in 2:ny)
   {
    y.grad <- pretty(min(x$linknodes[,i]):max(x$linknodes[,i]))
    y.grad[1] <- round(min(x$linknodes[,i]),2)
    y.grad[length(y.grad)] <- round(max(x$linknodes[,i]),2)
    lines(x=x$estimlink[,2*i],y=1:nsim,col=color[i])
    axis(ifelse(i%%2==0,4,2),at=loc.grad(y.grad,i),labels=y.grad,col=color[i], col.axis=color[i], cex.axis=0.8,line=(round(i/2)-1)*2)
   }
   legend(x=legend.loc, legend=x$Ynames, lty=1, col=color, bty="n", inset=c(0.05,0.05))
  }
  else
  {
   cat("Output can not be produced since the program stopped abnormally. \n")
  }
}

plot.linkfunction  <- function(x,legend.loc="topleft",...) UseMethod("plot.linkfunction")
