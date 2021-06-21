## ----setup, include = FALSE---------------------------------------------------
library(lcmm)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
  
)

## ---- , results='hide'--------------------------------------------------------
mhlme <- hlme(IST ~ I(age-age_init),random=~ I(age-age_init),subject="ID",data=paquid)
set.seed(1234)
mhlme2 <- hlme(IST ~ I(age-age_init),random=~ I(age-age_init),subject="ID",data=paquid,ng=2,
               mixture=~ I(age-age_init),classmb =~ CEP , B=random(mhlme))

## ---- comment=''--------------------------------------------------------------
summary(mhlme2)

## ---- , results='hide'--------------------------------------------------------
mhlme2perm <- permut(mhlme2, order=c(2,1))

## ---- comment=''--------------------------------------------------------------
summary(mhlme2)

## ---- comment=''--------------------------------------------------------------
xclass(mhlme2perm,mhlme2)

## ---- comment='', message = FALSE---------------------------------------------
predictClass(mhlme2, newdata=paquid[2:6,])

## ---- comment='',message = FALSE----------------------------------------------
predictRE(mhlme2, newdata=paquid[2:6,])

## ---- comment=''--------------------------------------------------------------
predss <- predictY(mhlme2, paquid[2:6,], marg=FALSE)
predss$pred

