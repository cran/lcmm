## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- results='hide',message=FALSE,warning=FALSE, echo=FALSE------------------
library(lcmm)
set.seed(1)

## ----comment=''---------------------------------------------------------------
library(NormPsy)
paquid$normMMSE <- normMMSE(paquid$MMSE)

## ----message=FALSE------------------------------------------------------------
paquid$age65 <- (paquid$age - 65)/10

## ---- comment='', results='hide'----------------------------------------------
m1 <- hlme(normMMSE ~ age65+I(age65^2)+CEP,random =~ age65+I(age65^2), subject = 'ID', data = paquid) # ng=1

## ----comment='', results='hide'-----------------------------------------------
#Estimation considering 2 classes : 
m2 <- hlme(normMMSE ~ age65+I(age65^2)+CEP, random =~ age65+I(age65^2), subject = 'ID', data = paquid, ng = 2, mixture=~age65+I(age65^2), B=m1)

## ---- comment='', results='hide'----------------------------------------------
m2b <- hlme(normMMSE ~ age65+I(age65^2)+CEP, random =~ age65+I(age65^2), subject = 'ID',data = paquid, ng = 2, mixture =~ age65+I(age65^2), B = c(0, 60, 40, 0, -4, 0, -10, 10, 212.869397, -216.421323,456.229910, 55.713775, -145.715516, 59.351000, 10.072221))

## ---- results='hide'----------------------------------------------------------
m2c <- hlme(normMMSE ~ age65+I(age65^2)+CEP, random =~ age65+I(age65^2),subject = 'ID', data = paquid, ng = 2, mixture =~ age65+I(age65^2),  B = random(m1))

## ---- results='hide', eval=FALSE----------------------------------------------
#  m2d <- gridsearch(hlme(normMMSE ~ age65+I(age65^2)+CEP,  random =~ age65+I(age65^2), subject = 'ID', data=paquid, ng = 2, mixture=~age65+I(age65^2)), rep=100, maxiter=30, minit=m1)

## ---- comment='', results='hide', eval=FALSE----------------------------------
#  m3g <- gridsearch(hlme(normMMSE ~ age65+I(age65^2)+CEP,  random =~ age65+I(age65^2), subject = 'ID', data=paquid, ng = 3, mixture=~age65+I(age65^2)), rep=100, maxiter=30, minit=m1)

## ---- comment='', results='hide', echo=FALSE----------------------------------
load("gridsearch_hlme.RData")

## ---- comment='', fig.height=3, fig.width=7-----------------------------------
summarytable(m1,m2,m2b,m2c, m2d  , m3g, which = c("G", "loglik", "conv", "npm", "AIC", "BIC", "SABIC", "entropy","ICL", "%class"))

summaryplot(m1, m2, m3g, which = c("BIC", "entropy","ICL"))

## ---- comment=''--------------------------------------------------------------
summary(m2d)

## ---- comment=''--------------------------------------------------------------
data_pred0 <- data.frame(age=seq(65,95,length.out=50),CEP=0)
data_pred1 <- data.frame(age=seq(65,95,length.out=50),CEP=1)
data_pred0$age65 <- (data_pred0$age - 65)/10
data_pred1$age65 <- (data_pred1$age - 65)/10

## ----comment=''---------------------------------------------------------------
pred0 <- predictY(m2d, data_pred0, var.time = "age")
pred1 <- predictY(m2d, data_pred1, var.time = "age")

## ----comment='', fig.height=4, fig.width=6------------------------------------
plot(pred0, col=c("red","navy"), lty=1,lwd=5,ylab="normMMSE",legend=NULL,  main="Predicted trajectories for normMMSE ",ylim=c(0,100))
plot(pred1, col=c("red","navy"), lty=2,lwd=3,legend=NULL,add=TRUE)
legend(x="topright",legend=c("class1 :","CEP-","CEP+","class2:","CEP-","CEP+"), col=c(rep("red",3),rep("navy",3)), lwd=2, lty=c(0,1,2,0,1,2), ncol=2, bty="n", cex = 0.7)

## ----comment='', fig.height=4, fig.width=6------------------------------------
predIC0 <- predictY(m2d, data_pred0, var.time = "age",draws=TRUE)
predIC1 <- predictY(m2d, data_pred1, var.time = "age",draws=TRUE)
plot(predIC0, col=c("deeppink","deepskyblue"), lty=1, lwd=2, ylab="normMMSE", main="Predicted trajectories for normMMSE", ylim=c(0,100), shades=TRUE)
plot(predIC1, col=c("deeppink","deepskyblue"), lty=2, lwd=2, ylab="normMMSE", main="Predicted trajectories for normMMSE", legend=NULL, ylim=c(0,100), shades=TRUE, add=TRUE)

## ----comment='',fig.height=4, fig.width=10------------------------------------
predG1 <- predictY(m1, data_pred0, var.time = "age")
predG3 <- predictY(m3g, data_pred0, var.time = "age")

par(mfrow=c(1,3))
plot(predG1, col=1, lty=1, lwd=2, ylab="normMMSE", legend=NULL, main="Predicted trajectories G=1",ylim=c(0,100))
plot(pred0, col=c("red","navy"), lty=1, lwd=2,ylab="normMMSE", legend=NULL, main="Predicted trajectories G=2", ylim=c(0,100))
plot(predG3, col=2:4, lty=1, lwd=2, ylab="normMMSE", legend=NULL, main="Predicted trajectories G=3", ylim=c(0,100))

## ---- eval=FALSE--------------------------------------------------------------
#  plot(m2d, cex.main=0.8)

## ---- comment='',fig.height=4, fig.width=6------------------------------------
plot(m2d, which="fit", var.time="age", marg=FALSE, shades = TRUE)

## ---- comment=''--------------------------------------------------------------
postprob(m2d) 

