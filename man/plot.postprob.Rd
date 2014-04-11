\name{plot.postprob}
\alias{plot.postprob.Jointlcmm}
\alias{plot.postprob}
\alias{plot.postprob.lcmm}
\alias{plot.postprob.hlme}
\alias{plot.postprob.multlcmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Histograms of the posterior class-membership probabilities
}
\description{
This function provides the histograms of the posterior class-membership probabilities stemmed from a \code{Jointlcmm}, \code{lcmm}, \code{hlme} or \code{multlcmm} object.
}
\usage{
\method{plot.postprob}{Jointlcmm}(x,\dots)
\method{plot.postprob}{lcmm}(x,\dots)
\method{plot.postprob}{hlme}(x,\dots)
\method{plot.postprob}{multlcmm}(x,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object inheriting from classes \code{lcmm}, \code{hlme}, \code{multlcmm} or \code{Jointlcmm}, representing respectively a fitted latent class linear mixed model, a more general latent class mixed model or a joint latent class model
}
  \item{\dots}{
other parameters to be passed through to plotting functions
}
}
\value{
Returns as many plots as the number of latent classes}
\author{
Cecile Proust-Lima, Viviane Philipps and Benoit Liquet
}
\note{
This function can only be used with latent class mixed models and joint latent class models that include at least 2 latent classes
}
\seealso{
\code{\link{Jointlcmm}}, \code{\link{lcmm}}, \code{\link{hlme}},\code{\link{postprob}}  
}
\examples{
data(data_hlme)
m<-lcmm(Y~Time*X1,mixture=~Time,random=~Time,classmb=~X2+X3,
subject='ID',ng=2,data=data_hlme,B=c(0.41,0.55,-0.18,-0.41,
-14.26,-0.34,1.33,13.51,24.65,2.98,1.18,26.26,0.97))  
plot.postprob(m) 
}
