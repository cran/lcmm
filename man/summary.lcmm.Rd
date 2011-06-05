\name{summary.lcmm}
\alias{summary.hlme}
\alias{summary.lcmm}
\alias{summary.Jointlcmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of a \code{hlme}, \code{lcmm} or \code{Jointlcmm} object
}
\description{
The function provides a summary of a \code{hlme}, \code{lcmm} or \code{Jointlcmm} estimation: model specification, goodness-of-fit statistics, maximum likelihood estimates and associated Wald test.
}
\usage{
\method{summary}{hlme}(object,\dots)
\method{summary}{lcmm}(object,\dots)
\method{summary}{Jointlcmm}(object,\dots)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object inheriting from classes \code{hlme}, \code{lcmm} or \code{Jointlcmm} representing a fitted latent class mixed-effects model or a Joint latent class mixed model.
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
\code{\link{hlme}},\code{\link{lcmm}},\code{\link{Jointlcmm}},
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{print}
