## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- results='hide',message=FALSE,warning=FALSE, echo=FALSE------------------
load("models_IRT.RData")

## ----message=F----------------------------------------------------------------
library(lcmm)
library(ggplot2)
library(ggpubr)
library(splines)
library(gridExtra)
library(dplyr)

## -----------------------------------------------------------------------------
str(simdataHADS)

## -----------------------------------------------------------------------------
demo <- simdataHADS %>% group_by(ID) %>% arrange(time) %>%
  filter(row_number()==1)
summary(demo)

## ----comment='', fig.height=4, fig.width=6------------------------------------
tempSim <- simdataHADS %>% group_by(ID) %>% arrange(time) %>% mutate(Visit=ifelse(time==time_entry,"Entry","Follow-up"))

p <- ggplot(tempSim, aes(x=time,fill=Visit,color=Visit)) + geom_histogram(binwidth=1,aes(y=..density..)) +
  labs(x="Months on the waiting list")
p + scale_color_grey(start = 0.1,
  end = 0.5)+scale_fill_grey(start = 0.1,
  end = 0.5) +
  theme_classic()

## -----------------------------------------------------------------------------
quantile(simdataHADS$time,probs=(0:10)/10)

## ----eval=F-------------------------------------------------------------------
#  modIRT_i <- multlcmm(hads_2 + hads_4 +hads_6 + hads_8 +hads_10+hads_12 + hads_14 ~ ns(time,knots=c(7,15),Boundary.knots = c(0,60))*grp,random=~1,data=simdataHADS,subject="ID",link="thresholds",methInteg="QMC",nMC=1000)
#  # use the estimates as initial values - the vector c(0,1,0,0,1,0,0,0,1) initializes the cholesky matrix of the random-effects
#  Binit <- c(modIRT_i$best[1:7],c(0,1,0,0,1,0,0,0,1),modIRT_i$best[8:length(modIRT_i$best)])
#  
#  modIRT <- multlcmm(hads_2 + hads_4 +hads_6 + hads_8 +hads_10+hads_12 + hads_14 ~ ns(time,knots=c(7,15),Boundary.knots = c(0,60))*grp,random=~1+ns(time,knots=c(7,15),Boundary.knots = c(0,60)),data=simdataHADS,subject="ID",link="thresholds",methInteg="QMC",nMC=1000, B=Binit)

## -----------------------------------------------------------------------------
summary(modIRT)

## ----comment='', fig.height=4, fig.width=6------------------------------------
datnew <- data.frame(time = seq(0,75,by=1))
datnew$grp <- 0
pIRT0 <- predictL(modIRT,datnew,var.time="time",confint = T)
datnew$grp <- 1
pIRT1 <- predictL(modIRT,datnew,var.time="time",confint = T)
plot(pIRT0,col=1,lwd=2,ylim=c(-1.5,1.5),legend=NULL,main="",ylab="latent depressive symptomatology",xlab="months since entry on the waiting list",type="l",bty="l",shades=T)
plot(pIRT1,add=T,col=4,lwd=2,shades=T)
legend(x="topleft",legend=c("dialysed","preemptive"),lty=c(1,1),col=c(1,4),lwd=2,bty="n")

## ----comment='', fig.height=4, fig.width=6------------------------------------
beta <- modIRT$best
t <- 0:72
Z <- cbind(rep(1,length(t)),ns(t,knots=c(7,15),Boundary.knots = c(0,60)))
chol <- matrix(0,ncol=4,nrow=4)
chol[upper.tri(chol, diag = T)] <- c(1,beta[7:15])
library(mvtnorm)
Lambda0 <- rmvnorm(10000,mean = Z%*%c(0,beta[1:3]),Z%*%t(chol)%*%chol%*%t(Z))
Lambda1 <- rmvnorm(10000,mean = Z%*%beta[4:7],Z%*%t(chol)%*%chol%*%t(Z))
Lambda <- data.frame(Lambda = as.vector(rbind(Lambda0,Lambda1)))
phist <- ggplot(Lambda,aes(x=Lambda))+ geom_density(color="grey", fill="grey") + theme_bw() +
  xlab("underlying depressive symptomatology") +xlim(-7,7)
phist

## -----------------------------------------------------------------------------
quantile(Lambda$Lambda,p=c(0.025,0.975))

## -----------------------------------------------------------------------------
## Parameters
nlevel <- 4
nitems <- 7
levels <- rep(nlevel,nitems)
npm <- length(modIRT$best)
seuils <- modIRT$best[(npm-(nlevel-1)*(nitems)+1):(npm)]
err <- abs(modIRT$best[(npm-(nlevel-1)*(nitems)-(nitems-1)):(npm-(nlevel-1)*(nitems))])
seuils
err
# Variance
Vseuils <- VarCov(modIRT)[(npm-(nlevel-1)*(nitems)+1):(npm),(npm-(nlevel-1)*(nitems)+1):(npm)]
Verr <- VarCov(modIRT)[(npm-(nlevel-1)*(nitems)-(nitems-1)):(npm-(nlevel-1)*(nitems)),(npm-(nlevel-1)*(nitems)-(nitems-1)):(npm-(nlevel-1)*(nitems))]

## -----------------------------------------------------------------------------
# generic function
location <- function(min,max,param,Vparam){
  loc <- param[1]
  se <- sqrt(Vparam[1,1])
  param2 <- rep(0,length(param))
  param2[1] <- 1
  if ((max-min)>1) {
    for (l in 1:(max-min-1)) {
      param2[l+1] <- 2*param[l+1]
      loc[l+1] <- loc[l] + param[1+l]^2
      se[l+1] <- sqrt(t(param2) %*%Vparam %*%param2)
    }
  }
  return(c(loc,se))
}
# application
ItemLoc <- sapply(1:nitems,function(k){location(min=0,max=nlevel-1,param=seuils[((nlevel-1)*(k-1)+1):((nlevel-1)*k)],Vparam=Vseuils[((nlevel-1)*(k-1)+1):((nlevel-1)*k),((nlevel-1)*(k-1)+1):((nlevel-1)*k)])})
colnames(ItemLoc) <- paste("Item",(1:nitems)*2)
ItemLoc <- ItemLoc[c(1,4,2,5,3,6),]
rownames(ItemLoc) <- rep(c("Threshold","SE"),nlevel-1)
discrimination <- 1/abs(err)
sediscr <- diag(err^(-2))%*%Verr%*%diag(err^(-2))

## -----------------------------------------------------------------------------
t(rbind(ItemLoc,discrimination,Se=sqrt(diag(sediscr))))

## ----comment='', fig.height=4, fig.width=7------------------------------------
## computations
info_modIRT <- ItemInfo(modIRT, lprocess=seq(-6,6,0.1))

meaning <- c("Enjoy","Laugh","Cheerful" ,"Slow" ,"Appearance" ,"Looking Forward" ,"Leisure")
items <- paste("hads", seq(2,14,2), sep="_")

## automatic graph
par(mfrow=c(2,4), mar=c(3,2,2,1), mgp=c(2,0.5,0))
for(k in 1:7)
{     
 plot(info_modIRT, which="LevelProb", outcome=items[k],
      main=paste("Item",2*k,"-",meaning[k]), lwd=2, legend=NULL)
}
plot(0,axes=FALSE, xlab="", ylab="", col="white")
legend("center", legend=paste("Level",0:3), col=1:4, lwd=2, box.lty=0)

## ----comment='', fig.height=4, fig.width=7------------------------------------
## graph with ggplot
p <- NULL
for (k in 1:7){
ICC  <- info_modIRT$LevelInfo[which(info_modIRT$LevelInfo[,1]==items[k]),]
p[[k]] <- (ggplot(ICC)
      + geom_line(aes(x = Lprocess, y = Prob, group = Level,color=Level), show.legend = F,alpha = 1,size=1.2)
      # + stat_smooth(aes(x = time, y = hads_scorea), method = "loess", size = 0.75)
      + theme_bw()
      + labs(title=paste("Item",2*k,"-",meaning[k]))
      + xlab("construct")
      + ylim(0,1)
      + ylab("Probability of item level"))
}
p[[8]] <- (ggplot(ICC)
      + geom_line(aes(x = Lprocess, y = Prob, group = Level,color=Level),alpha = 1,size=1.2)
      + theme_bw()
)
legend <- get_legend(p[[8]])
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],as_ggplot(legend),ncol=4)

## ----comment='', fig.height=4, fig.width=7------------------------------------
expe <- predictYcond(modIRT,lprocess = seq(-6,6,by=0.1))
# via the internal plot function 
plot(expe, xlab="underlying depressive symptomatology", main="Item Expectation Curves",
     legend=paste("Item",(1:nitems)*2), lwd=2)

## ----comment='', fig.height=4, fig.width=6------------------------------------
# via ggplot
j <- table(expe$pred$Yname)[1]
expe$pred$item <- as.factor(c(rep(2,j),rep(4,j),rep(6,j),rep(8,j),rep(10,j),rep(12,j),rep(14,j)))
p <- (ggplot(expe$pred)
      + geom_line(aes(x = Lprocess, y = Ypred, group=item,color=item),alpha = 1,size=1.2)
      + theme_bw()
      + xlab("underlying depressive symptomatology")
      + ylim(0,3)
      + ylab("Item Expectation"))
p

## ----message=F,warning=F,comment='', fig.height=4, fig.width=7----------------
par(mfrow=c(2,4), mar=c(3,2,2,1), mgp=c(2,0.5,0))
for(k in 1:7)
{     
 plot(info_modIRT, which="LevelInfo", outcome=items[k],
 main=paste("Item",2*k,"-",meaning[k]), lwd=2, legend=NULL, ylim=c(0,1.3))
}
plot(0,axes=FALSE, xlab="", ylab="", col="white")
legend("center", legend=paste("Level",0:3), col=1:4, lwd=2, box.lty=0)

## ----message=F,warning=F,comment='', fig.height=4, fig.width=6----------------
plot(info_modIRT, which="ItemInfo", lwd=2, legend.loc="topleft")

## ----eval=F-------------------------------------------------------------------
#  head(datnew)
#  datnew$grp <- 0
#  ns0 <- predictY(modIRT,var.time = "time",newdata=datnew,methInteg = 1,nsim=2000,draws=T)
#  datnew$grp <- 1
#  ns1 <- predictY(modIRT,var.time = "time",newdata=datnew,methInteg = 1,nsim=2000,draws=T)

## ----comment='', fig.height=4, fig.width=6------------------------------------
par(mfrow=c(2,4), mar=c(3,2,2,1), mgp=c(2,0.5,0))
for(k in 1:7){
plot(ns0,outcome = k,shades = T,ylim=c(0,3),bty="l",legend=NULL,main=paste("Item",2*k,"-",meaning[k]),ylab="Item level",xlab="months on the waiting list")
plot(ns1,outcome=k,shades=T,add=T,col=2)
}

## ----eval=F-------------------------------------------------------------------
#  # initialization of the parameter vector for faster convergence
#  npm <- length(modIRT$best)
#  Binit <- c(modIRT$best[1:7],rep(0,(nitems-1)),modIRT$best[(npm-nlevel*nitems-9+1):npm])
#  # estimation
#  modIRT_DIFg <- multlcmm(hads_2 + hads_4 +hads_6 + hads_8 +hads_10+hads_12 + hads_14 ~ ns(time,knots=c(7,15),Boundary.knots = c(0,60))*(grp) +contrast(grp),random=~1+ns(time,knots=c(7,15),Boundary.knots = c(0,60)),data=simdataHADS,subject="ID",link="thresholds",methInteg="QMC",nMC=1000,B=Binit)

## -----------------------------------------------------------------------------
sumDIF <- summary(modIRT_DIFg)
sumDIF[,2]

## -----------------------------------------------------------------------------
WaldMult(modIRT_DIFg,pos=c(8:13))

## -----------------------------------------------------------------------------
sum <- summary(modIRT_DIFg)[10,]
c(sum[1],sum[1]- qnorm(0.975)*sum[2],sum[1]+ qnorm(0.975)*sum[2])

## ----eval=F-------------------------------------------------------------------
#  # initialization of the parameter vector for faster convergence
#  npm <- length(modIRT$best)
#  Binit <- c(modIRT$best[1:7],rep(0,3*(nitems-1)),modIRT$best[(npm-nlevel*nitems-9+1):npm])
#  # estimation
#  modIRT_DIFt <- multlcmm(hads_2 + hads_4 +hads_6 + hads_8 +hads_10+hads_12 + hads_14 ~ ns(time,knots=c(7,15),Boundary.knots = c(0,60))*(grp) + contrast(ns(time,knots=c(7,15),Boundary.knots = c(0,60))),random=~1+ns(time,knots=c(7,15),Boundary.knots = c(0,60)),data=simdataHADS,subject="ID",link="thresholds",methInteg="QMC",nMC=1000,B=Binit)

## -----------------------------------------------------------------------------
summary(modIRT_DIFt)

## -----------------------------------------------------------------------------
## Pvalue for the last contrast
b <- coef(modIRT_DIFt)
v <- vcov(modIRT_DIFt)
A <- rbind(c(rep(0,7), rep(-1,6), rep(0,49)),
	   c(rep(0,7+6), rep(-1,6), rep(0,49-6)),
	   c(rep(0,7+12), rep(-1,6), rep(0,49-12)))
w <- t(A%*%b) %*% solve(A%*%v%*%t(A)) %*% A%*%b
DIF14 <- 1-pchisq(w, df=nrow(A)) # p=0.3722833

## -----------------------------------------------------------------------------
## pvalues for all the items including the last one
DIF <- cbind(seq(2,14,by=2),c(sapply(1:6,function(k){WaldMult(modIRT_DIFt,pos=c(7+k,(7+6+k),(7+2*6+k)))[2]}),DIF14))
colnames(DIF) <- c("item","pvalue")
DIF   

## ----eval=F-------------------------------------------------------------------
#  head(datnew)
#  datnew$grp <- 0
#  ns0DIFt <- predictY(modIRT_DIFt,var.time = "time",newdata=datnew,methInteg = 1,nsim=2000,draws=T)
#  datnew$grp <- 1
#  ns1DIFt <- predictY(modIRT_DIFt,var.time = "time",newdata=datnew,methInteg = 1,nsim=2000,draws=T)

## ----comment='', fig.height=6, fig.width=6------------------------------------
par(mfrow=c(3,3), mar=c(3,2,2,1), mgp=c(2,0.5,0))
for(k in 1:7){
plot(ns0,outcome = k,shades = T,ylim=c(0,3),bty="l",legend=NULL,main=paste("Item",2*k,"-",meaning[k]),ylab="Item level",xlab="months on the waiting list",color=1,lwd=2,xlim=c(0,50))
plot(ns0DIFt,outcome=k,shades=T,lty=2,add=T,col=1,lwd=2)
plot(ns1,outcome=k,shades=T,add=T,col=4,lwd=2)
plot(ns1DIFt,outcome=k,shades=T,add=T,col=4,lty=2,lwd=2)
legend("top",legend=paste("(RS overall test: p = ",round(DIF[k,2],digits = 3),")",sep=""),bty="n")
}
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', main ='')
legend(x="top",legend=c("dialysed","pre-emptive"),lty=c(1,1),col=c(1,4),lwd=2,bty="n")
legend(x="bottom",legend=c("without RS","with RS"),lty=c(1,2),col=c("gray","gray"),lwd=2,bty="n")

