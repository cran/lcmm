## ----setup, include = FALSE---------------------------------------------------
library(lcmm)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
  
)

## ---- eval=FALSE--------------------------------------------------------------
#  hlme(fixed, mixture, random, subject, classmb, ng = 1, idiag = FALSE, nwg = FALSE, cor = NULL, data, B, convB = 0.0001, convL = 0.0001, convG = 0.0001, prior, maxiter = 500, subset = NULL, na.action = 1, posfix = NULL)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("lcmm")

## ----results='hide',message=FALSE,warning=FALSE-------------------------------
library(lcmm)

## ----comment=''---------------------------------------------------------------
head(paquid)

## ----comment=''---------------------------------------------------------------
summary(paquid)

## ---- results='hide',message=FALSE,warning=FALSE, echo=FALSE------------------
library(lcmm)

## ----comment='', results='hide'-----------------------------------------------
library(NormPsy)

## ----comment='', fig.height=4, fig.width=6------------------------------------
paquid$normMMSE <- normMMSE(paquid$MMSE)
par(mfrow=c(1,2))
hist(paquid$MMSE, cex.main=0.8, cex.lab=0.8)
hist(paquid$normMMSE, cex.main=0.8, cex.lab=0.8)

## ----comment='', fig.height=4, fig.width=6------------------------------------
library(lattice)
color <- paquid$ID
xyplot(normMMSE ~ age, paquid, groups = ID, col=color, lwd=2, type="l")

