## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>")

## ---- include=FALSE-----------------------------------------------------------
library(lcmm)


## ----message=FALSE------------------------------------------------------------
paquid$time <- (paquid$age - paquid$age_init)/10
paquid$age75 <- (paquid$age_init - 75)/10

## ---- comment='', results='hide', echo=FALSE----------------------------------
load("models_mult.RData")

## ----message=FALSE, results='hide',eval=FALSE---------------------------------
#  mult_lin <- multlcmm(MMSE + IST + BVRT ~ age75 + male + time + I(time^2/10) + contrast(male), random =~ time +I(time^2 / 10), subject='ID', data = paquid, randomY = TRUE, cor = BM(time))

## ----hist,message=FALSE,fig.height=6,fig.width=7------------------------------
par(mfrow=c(2,2))
hist(paquid$MMSE, breaks=31,main="MMSE distribution")
hist(paquid$IST, breaks=41,main="IST distribution")
hist(paquid$BVRT, breaks=16, main="BVRT distribution")

## ----test2,message=FALSE, results='hide',eval=FALSE---------------------------
#  # Example with Beta
#  mult_beta <- multlcmm(MMSE + IST + BVRT ~ age75 + male + time + I(time^2/10) + contrast(male), random =~ time + I(time^2/10), subject='ID', data = paquid, randomY = TRUE, cor = BM(time), link = 'beta')

## ----test73,message=FALSE, results='hide',eval=FALSE--------------------------
#  # different number of knots in splines
#  mult_betaspl <- multlcmm(MMSE + IST + BVRT ~ age75 + male + time + I(time^2/10) + contrast(male), random =~time + I(time^2/10), subject='ID', data = paquid, randomY = TRUE, cor = BM(time), link = c('beta','3-quant-splines','3-quant-splines'))

## ----test71,message=FALSE, results='hide',eval=FALSE--------------------------
#  # with splines
#  mult_splines <- multlcmm(MMSE + IST + BVRT ~ age75 + male + time + I(time^2/10) + contrast(male), random =~time + I(time^2/10), subject='ID', data = paquid, randomY = TRUE, cor = BM(time), maxiter=50, link = c('3-quant-splines'))

## ----test27,message=FALSE-----------------------------------------------------
mult_splines$best

## ----test72,message=FALSE, results='hide',eval=FALSE--------------------------
#  # with splines
#  mult_splines2 <- multlcmm(MMSE + IST + BVRT ~ age75 + male + time + I(time^2/10) + contrast(male), random =~time + I(time^2/10), subject='ID', data = paquid, randomY = TRUE, cor = BM(time), maxiter=50, link = c('3-quant-splines'), posfix=21, B=mult_splines$best)

## ----test5,message=FALSE,comment=''-------------------------------------------
summarytable(mult_lin,mult_beta,mult_betaspl,mult_splines2, which =c("loglik", "conv", "npm", "AIC"))

## ----testprint2,message=FALSE, fig.height=3, fig.width=5----------------------
par(mfrow=c(1,1))
col <- rainbow(4)
plot(mult_splines2, which = "linkfunction", col = c(col[2],col[3],col[4]), lwd =1,lty=4)
plot(mult_lin,which="linkfunction", col = c(col[2],col[3],col[4]), lwd = 1,lty=2,add=TRUE)
plot(mult_beta,which="linkfunction", col = c(col[2],col[3],col[4]), lwd = 2,lty=3,add=TRUE)
plot(mult_betaspl,which="linkfunction", col = c(col[2],col[3],col[4]), lwd = 1,lty=1,add=TRUE)
legend(x="bottomright",lty=c(2,3,4,1),legend=c("linear","beta","splines","beta/splines"),bty="n")

## ----testprint,message=FALSE, fig.height=3, fig.width=5-----------------------
col <- rainbow(4)
plot(mult_betaspl, which = "linkfunction", col = c(col[2],col[3],col[4]), lwd = 2)

## ----message=FALSE,echo=FALSE,  fig.height=3, fig.width=5---------------------
col <- rainbow(5)
CIlink <- predictlink(mult_betaspl)
plot(CIlink, col = c(col[2],col[3],col[4]), lwd = 2, shades = TRUE)

## ----message=FALSE, comment=''------------------------------------------------
summary(mult_betaspl) 

## ---- comment=''--------------------------------------------------------------
VarExpl(mult_betaspl, data.frame(time=0))

## ----message=FALSE,eval=TRUE,  fig.height=6, fig.width=7----------------------
datnew <- data.frame(time=seq(0.08,2.2,length=100))
datnew$age_init<-seq(65,95, length=100)
datnew$age75 <- ((datnew$age_init - 75)/10)
datnew$male <- 0
predict_women<-predictY(mult_betaspl,newdata=datnew,var.time='time',draws=TRUE)

datnew$male <- 1
predict_men <- predictY(mult_betaspl,newdata=datnew,var.time='time',draws=TRUE) 
par(mfrow=c(2,2))
plot(predict_women, lwd=c(2,1), type='l', col=6, ylim=c(0,30), xlab='time since entry (in decades)', ylab='Marker', bty='l', legend=NULL, shades=TRUE, outcome = 1, main='Mean predicted trajectory for MMSE')
plot(predict_men, lwd=c(2,1), type='l', col=3, shades=TRUE, outcome = 1, add=TRUE)
legend(1.5, 30, legend=c("Women", "Men"), col=c(6,3), lty=1:2, cex=0.8,bty="n")
plot(predict_women, lwd=c(2,1), type='l', col=6, ylim=c(0,40), xlab='time since entry (in decades)', ylab='Marker', bty='l', legend=NULL, shades=TRUE, outcome = 2, main='Mean predicted trajectory for IST') 
plot(predict_men, lwd=c(2,1), type='l', col=3, shades=TRUE, outcome = 2, add=TRUE)
plot(predict_women, lwd=c(2,1), type='l', col=6, ylim=c(0,15), xlab='time since entry (in decades)', ylab='Marker', bty='l', legend=NULL, shades=TRUE, outcome = 3, main='Mean predicted trajectory for BVRT') 
plot(predict_men, lwd=c(2,1), type='l', col=3, shades=TRUE, outcome = 3, add=TRUE)

## ----message=FALSE,eval=TRUE,  fig.height=6, fig.width=7----------------------
plot(mult_betaspl, cex.main=0.8)

## ----message=FALSE,eval=TRUE,  fig.height=6, fig.width=7----------------------
par(mfrow=c(2,2))
plot(mult_betaspl, which="fit", var.time="time", bty="l", xlab="time since entry (in decades)", cex.lab=1.1, break.times=8, ylab="latent process", lwd=2, marg=FALSE, ylim=c(-2,0.0), xlim=c(0.1,2), shades = TRUE, outcome = 1, col=2, main="MMSE predictions vs observations") 

plot(mult_betaspl, which="fit", var.time="time", bty="l", xlab="time since entry (in decades)", cex.lab=1.1, break.times=8, ylab="latent process", lwd=2, marg=FALSE, ylim=c(-2,0.3), xlim=c(0.1,2), shades = TRUE, outcome = 2, col=3, main="IST predictions vs observations")

plot(mult_betaspl, which="fit", var.time="time", bty="l", xlab="time since entry (in decades)", cex.lab=1.1, break.times=8, ylab="latent process", lwd=2, marg=FALSE, ylim=c(-1.5,0.5), xlim=c(0.1,2), shades = TRUE, outcome = 3, col=4, main="BVRT predictions vs observations") 


