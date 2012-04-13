\name{plot.postprob}
\alias{plot.postprob.Jointlcmm}
\alias{plot.postprob}
\alias{plot.postprob.lcmm}
\alias{plot.postprob.hlme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Histograms of the posterior class-membership probabilities
}
\description{
This function provides the histograms of the posterior class-membership probabilities stemmed from a \code{Jointlcmm} object.
}
\usage{
\method{plot.postprob}{Jointlcmm}(x,\dots)
\method{plot.postprob}{lcmm}(x,\dots)
\method{plot.postprob}{hlme}(x,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object inheriting from classes \code{lcmm}, \code{hlme} or \code{Jointlcmm}, representing respectively a fitted latent class linear mixed model, a more general latent class mixed model or a joint latent class model
}
  \item{\dots}{
other parameters to be passed through to plotting functions
}
}
\value{
Returns as many plots as the number of latent classes}
\author{
Cecile Proust-Lima and Benoit Liquet
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
subject='ID',ng=2,data=data_hlme)  
plot.postprob(m) 
}
