\name{postprob.hlme}
\alias{postprob.hlme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior classification stemmed from a LCMM estimation
}
\description{
This function provides informations about the posterior classification stemmed from a \code{hlme} object.
}
\usage{
\method{postprob}{hlme}(x,\dots)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{x}{an object inheriting from class \code{hlme}, representing a fitted latent class
          linear mixed-effects model}
 \item{\dots}{
further arguments to be passed to or from other methods.
          They are ignored in this function.
}
}

\details{
This function provides the number of subjects classified a posteriori in each latent class, and the classification table that contains the mean of the posterior probability of belonging to each latent class over the subjects classified in each of the latent classes. This table aims at evaluating the quality of the posterior classification
}
\author{
Cecile Proust-Lima and Benoit Liquet
}
\note{
This function can only be used with latent class mixed models that include at least 2 latent classes
}
\seealso{
\code{\link{plot.postprob.hlme}}, \code{\link{hlme}}, \code{\link{plot.predict.hlme}}
}
\examples{
data(data_hlme)
m<-hlme(Y~Time+X1+X1_time,mixture=~Time,random=~Time,classmb=~X2+X3,
 subject=ID,ng=2,data=data_hlme,B=c(0,0,0,30,25,0,-1,0,0,5,0,1,1))
postprob.hlme(m)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

