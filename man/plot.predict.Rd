\name{plot.predict}
\alias{plot.predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Marginal class-specific predictions plot
}
\description{
This function is generic; method functions can be written to
     handle specific classes of objects. For the moment, the only class which already has
     methods for this function is \code{hlme}.
}
\usage{
plot.predict(x,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object inheriting from class \code{hlme}, representing a fitted latent class
          linear mixed-effects model
}
  \item{\dots}{
further arguments to be passed to or from other methods.
          They are ignored in this function.
}
}
\value{
Returns a plot}
\author{
Cecile Proust-Lima and Benoit Liquet
}
\seealso{
\code{\link{hlme}}, \code{\link{plot.postprob}}, \code{\link{postprob}}, \code{\link{plot.predict.hlme}} 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
