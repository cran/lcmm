\name{plot.postprob}
\alias{plot.postprob}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Histograms of the posterior class-membership probabilities
}
\description{
This function is generic; method functions can be written to
     handle specific classes of objects. For the moment, the only class which already has
     methods for this function is \code{hlme}.
}
\usage{
plot.postprob(x,\dots)
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
\code{\link{plot.postprob.hlme}}, \code{\link{hlme}}, \code{\link{postprob}}, \code{\link{plot.predict}}  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
