## ----setup, include = FALSE---------------------------------------------------
library(lcmm)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- comment=''--------------------------------------------------------------
summary(paquid$CESD)

## ---- fig.height=4, fig.width=6-----------------------------------------------
hist(paquid$CESD, breaks=50)

## ---- comment=''--------------------------------------------------------------
#We recenter and scale the time variable "age" in order to avoid numerical problems
paquid$age65 <- (paquid$age-65)/10

## ---- results='hide'----------------------------------------------------------
mpreH <- lcmm(CESD ~ age65 + I(age65^2), random = ~ age65 + I(age65^2), subject = 'ID', data=paquid, link = '5-quant-splines') 

## ---- comment=''--------------------------------------------------------------
head(mpreH$pred)

## ---- comment=''--------------------------------------------------------------
paquid$normCESD <- NULL 
paquid$normCESD[!is.na(paquid$CESD)] <- mpreH$pred$obs

## ---- comment=''--------------------------------------------------------------
summary(paquid[,c("CESD","normCESD")])

## ---- fig.height=4, fig.width=7, comment=''-----------------------------------
par(mfrow=c(1,2))
hist(paquid$CESD, breaks=50, cex.main=0.9, main="Distribution of CESD")
hist(paquid$normCESD, breaks=50, cex.main=0.9, main="Distribution of normCESD") 

## ---- results='hide'----------------------------------------------------------
normCESD <- hlme(normCESD ~ age65*male, random = ~ age65, subject = 'ID', data=paquid)

## ---- eval=FALSE, comment=''--------------------------------------------------
#  plot(normCESD, cex.main=0.8)

## ---- results='hide'----------------------------------------------------------
CESD <- hlme(CESD ~ age65*male, random = ~ age65, subject = 'ID', data=paquid)

## ---- eval=FALSE, comment=''--------------------------------------------------
#  plot(CESD, cex.main=0.8)

## ---- comment=''--------------------------------------------------------------
m <- mean(paquid$normCESD[(paquid$visit==0) & (!is.na(paquid$normCESD))])
s <- sd(paquid$normCESD[(paquid$visit==0) & (!is.na(paquid$normCESD))])
paquid$ZnormCESD <- (paquid$normCESD - m)/s

## ---- comment=''--------------------------------------------------------------
min <- min(paquid$normCESD[!is.na(paquid$normCESD)])
max <- max(paquid$normCESD[!is.na(paquid$normCESD)])
paquid$normCESD100 <- (paquid$normCESD - min)/(max-min)*100
summary(paquid$normCESD100)

## ---- results='hide', comment=''----------------------------------------------
m1 <- hlme(normCESD100 ~ age65*male + CEP*age65 + age_init, random=~age65, subject='ID',data=paquid)
summary(m1)

## ---- results='hide', comment=''----------------------------------------------
paquid$time <- paquid$age - paquid$age_init
m2 <- hlme(normCESD100 ~ time*male + CEP*time + age_init, random=~time, subject='ID', data=paquid)
summary(m2)

