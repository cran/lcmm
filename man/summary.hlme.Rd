\name{summary.hlme}
\alias{summary.hlme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of a hlme object
}
\description{
The function provides a summary of the \code{hlme} estimation: model specification, goodness-of-fit statistics, maximum likelihood estimates and associated Wald test.
}
\usage{
\method{summary}{hlme}(object,\dots)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object inheriting from class \code{hlme}, representing a fitted latent class
          linear mixed-effects model
}
  \item{\dots}{
further arguments to be passed to or from other methods.
          They are ignored in this function.
}
}
\author{
Cecile Proust-Lima and Benoit Liquet
}
\seealso{
\code{\link{hlme}},\code{\link{postprob.hlme}},\code{\link{plot.postprob.hlme}},
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{print}
