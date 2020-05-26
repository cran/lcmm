## ----setup, include = FALSE---------------------------------------------------
library(lcmm)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
  
)

## ---- include=FALSE,echo=FALSE------------------------------------------------
library(lcmm)
paquid$age65 <- (paquid$age - 65)/10

## ----message=FALSE, results='hide'--------------------------------------------
mlin <- lcmm(CESD ~ age65*male, random=~ age65, subject='ID', data=paquid) #link= linear

## ----message=FALSE, results='hide'--------------------------------------------
mlin2 <- hlme(CESD ~ age65*male, random=~ age65, subject='ID', data=paquid) #link= linear

## ---- comment=''--------------------------------------------------------------
mlin$loglik
mlin2$loglik
mlin$best
mlin2$best

## ----test2,message=FALSE, results='hide'--------------------------------------
mbeta <- lcmm(CESD ~ age65*male, random=~ age65, subject='ID', data=paquid, link='beta')

## ----test7,message=FALSE, results='hide'--------------------------------------
mspl <- lcmm(CESD ~ age65*male, random=~ age65, subject='ID', data=paquid, link='splines')

## ----test8,message=FALSE, results='hide'--------------------------------------
mspl5q <- lcmm(CESD ~ age65*male, random=~ age65, subject='ID', data=paquid, link='5-quant-splines')

## ----test5,message=FALSE,comment=''-------------------------------------------
summarytable(mlin,mbeta,mspl,mspl5q,which = c("loglik", "conv", "npm", "AIC"))

## ----testprint,message=FALSE, fig.height=4, fig.width=6-----------------------
col <- rainbow(5)
plot(mlin, which="linkfunction", bty='l', ylab="CES-D", col=col[1], lwd=2, xlab="underlying latent process")
plot(mbeta, which="linkfunction", add=TRUE, col=col[2], lwd=2)
plot(mspl, which="linkfunction", add=TRUE, col=col[3], lwd=2)
plot(mspl5q, which="linkfunction", add=TRUE, col=col[4], lwd=2)
legend(x="topleft", legend=c("linear", "beta","splines (5equidistant)","splines (5 at quantiles)"), lty=1, col=col, bty="n", lwd=2)

## ----message=FALSE,  fig.height=4, fig.width=6--------------------------------
linkspl5q <- predictlink(mspl5q,ndraws=2000)
plot(linkspl5q, col=col[4], lty=2, shades=TRUE)
legend(x="left", legend=c("95% confidence bands","for splines at quantiles"),lty=c(2,NA), col=c(col[4],NA), bty="n", lwd=1, cex=0.8)

## ---- results='hide', eval=FALSE----------------------------------------------
#  mthresholds <- lcmm(HIER ~ age65*male, random=~ age65, subject='ID', data=paquid, link='thresholds')

## ----message=FALSE, comment=''------------------------------------------------
summary(mspl5q) 

## ----message=FALSE,eval=TRUE--------------------------------------------------
datnew <- data.frame(age=seq(65,95,length=100))
datnew$age65 <- (datnew$age - 65)/10
datnew$male <- 0
women <- predictY(mspl5q, newdata=datnew, var.time="age", draws=TRUE)
datnew$male <- 1
men <- predictY(mspl5q, newdata=datnew, var.time="age", draws=TRUE)

## ----message=FALSE,  fig.height=4, fig.width=6--------------------------------
plot(women, lwd=c(2,1), type="l", col=6, ylim=c(0,20), xlab="age in year",ylab="CES-D",bty="l", legend=NULL, shades = TRUE)
plot(men, add=TRUE, col=4, lwd=c(2,1), shades=TRUE)
legend(x="topleft", bty="n", ncol=2, lty=c(1,1,2,2), col=c(6,4,6,4), legend=c("women","men", "95% CI", "95% CI"), lwd=c(2,2,1,1)) 

## ----message=FALSE,eval=TRUE,  fig.height=5, fig.width=8----------------------
plot(mspl5q, cex.main=0.9)

## ----message=FALSE,eval=TRUE,  fig.height=4, fig.width=6----------------------
plot(mspl5q, which="fit", var.time="age65", bty="l", xlab="(age-65)/10", break.times=8, ylab="latent process", lwd=2, marg=FALSE, ylim=c(-1,2), shades=TRUE, col=2)

