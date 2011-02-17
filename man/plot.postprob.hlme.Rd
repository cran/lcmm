\name{plot.postprob.hlme}
\alias{plot.postprob.hlme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Histograms of the posterior class-membership probabilities
}
\description{
This function provides the histograms of the posterior class-membership probabilities stemmed from a \code{hlme} object.
}
\usage{
\method{plot.postprob}{hlme}(x,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object inheriting from class \code{hlme}, representing a fitted latent class
          linear mixed-effects model
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
This function can only be used with latent class mixed models that include at least 2 latent classes
}
\seealso{
\code{\link{postprob.hlme}}, \code{\link{hlme}}, \code{\link{plot.predict.hlme}}  
}
\examples{
data(data_hlme)
m<-hlme(Y~Time+X1+X1_time,mixture=~Time,random=~Time,classmb=~X2+X3,subject=ID,ng=2,data=data_hlme,B=c(0,0,0,30,25,0,-1,0,0,5,0,1,1))

plot.postprob.hlme(m)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
