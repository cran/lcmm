\name{lcmm-package}
\alias{lcmm-package}
\alias{lcmm}
\docType{package}
\title{
Estimation of Latent Class Mixed Models and Joint Latent Class Mixed Models for quantitative, bounded quantitative, discrete and ordinal outcomes
}
\description{
This package provides functions for the estimation of latent class mixed models (LCMM) and joint latent class mixed models using a maximum likelihood method.
}
\details{
\tabular{ll}{
Package: \tab lcmm \cr
Type: \tab Package \cr
Version: \tab 1.2 \cr
Date: \tab 2011-02-17 \cr
License: \tab GPL (>=2.0) \cr
LazyLoad: \tab yes \cr
}
The package includes for the moment only the estimation of latent class linear mixed models for Gaussian longitudinal outcomes using \code{hlme} function. 
}
\author{
Cecile Proust-Lima, Amadou Diakite and Benoit Liquet

Maintainer: <cecile.proust@isped.u-bordeaux2.fr>
}
\references{
Verbeke G and Lesaffre E (1996). A linear mixed-effects model with heterogeneity in
the random-effects population. Journal of the American Statistical Association 91, 217-21

Muthen B and Shedden K (1999). Finite mixture modeling with mixture outcomes using the EM algorithm. Biometrics 55, 463-9

Proust C and Jacqmin-Gadda H. Estimation of linear mixed models with a mixture of distribution for the random-effects. Comput Methods Programs Biomed 78:165-73
}
\examples{
data(data_hlme)
m<-hlme(Y~Time+X1+X1_time,mixture=~Time,random=~Time,classmb=~X2+X3,subject=ID,ng=2,data=data_hlme,B=c(0,0,0,30,25,0,-1,0,0,5,0,1,1))
summary(m)
plot(m)
##posterior classification
postprob(m)
plot.postprob(m)
##class-specific predicted trajectories
newdata<-data.frame(intercept=rep(1,100),Time=seq(0,5,length=100),X1=rep(0,100),X1_time=rep(0,100),X2=rep(0,100),X3=rep(0,100))
plot.predict.hlme(m,newdata,"Time","right")
}
\keyword{package}

