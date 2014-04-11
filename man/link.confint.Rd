\name{link.confint}
\alias{link.confint.lcmm}
\alias{link.confint.multlcmm}
\alias{link.confint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confidence intervals for the estimated link functions from \code{lcmm} and \code{multlcmm}
}
\description{
This function provides 95\% confidence intervals around the estimated transformation given in estimlink attribute of \code{lcmm} and \code{multlcmm} objects.
It can also be used to evaluate the link functions at other values than those given in attribute \code{estimlink} of \code{lcmm} or \code{multlcmm} object.
Furthermore, a graph of the link functions, representing the relation between the underlying latent process and the observed outcomes, and the confidence intervals can be plotted.
}
\usage{
\method{link.confint}{lcmm}(x,ndraws=2000,Yvalues,plot=TRUE,...)
\method{link.confint}{multlcmm}(x,ndraws=2000,Yvalues,plot=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object inheriting from classes \code{lcmm} or \code{multlcmm}.
}
  \item{ndraws}{
the number of draws that should be generated to approximate the posterior distribution of the transformed values. By default, ndraws=2000.
}
  \item{Yvalues}{
a vector (for a \code{lcmm} object) or a matrix (for a \code{multlcmm} object) containing the values at which to compute the transformation(s).
Default to the values in \code{x$estimlink}.
}
  \item{plot}{
logical indicating if the link functions and confidence intervals should be plotted. Default to TRUE.
}
  \item{\dots}{
other parameters to be passed through to plotting functions
}
}
\value{
For a \code{lcmm} object, returns a data frame containing the values at which the transformation is evaluated,
the transformed values and the lower and the upper limits of the confidence intervals.

For a \code{multlcmm} object, returns a data frame containing the indicator of the outcome, the values at which the transformations are evaluated,
the transformed values and the lower and the upper limits of the confidence intervals.

}
\author{
Cecile Proust-Lima and Viviane Philipps
}
\seealso{
\code{\link{lcmm}}, \code{\link{multlcmm}}, \code{\link{plot.linkfunction}}
}
\examples{
 \dontrun{
data(data_Jointlcmm)

## Univariate mixed model with splines link funciton
m14<-lcmm(Ydep2~Time+I(Time^2),random=~Time,subject='ID',ng=1,
data=data_Jointlcmm,link="5-manual-splines",intnodes=c(10,20,25))
##Transformed values of several scores and their confidence intervals
transf.m14 <- link.confint(m14,ndraws=2000,Yvalues=c(0,1,7:30),plot=TRUE)

## Multivariate mixed model with splines link functions
m1 <- multlcmm(Ydep1+Ydep2~1+Time*X2+contrast(X2),random=~1+Time,
subject="ID",randomY=TRUE,link=c("4-manual-splines","3-manual-splines"),
intnodes=c(8,12,25),data=data_Jointlcmm,
B=c(-1.071, -0.192,  0.106, -0.005, -0.193,  1.012,  0.870,  0.881,
  0.000,  0.000, -7.520,  1.401,  1.607 , 1.908,  1.431,  1.082,
 -7.528,  1.135 , 1.454 , 2.328, 1.052))
##Confidence intervals for the transformed values (given in m1$estimlink)
transf.m1 <- link.confint(m1,ndraws=2000,plot=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
