\name{print.lcmm}
\alias{print.hlme}
\alias{print.lcmm}
\alias{print.Jointlcmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Short summary of a \code{hlme}, \code{lcmm} or \code{Jointlcmm} object
}
\description{
The function provides a summary of a \code{hlme}, \code{lcmm} or \code{Jointlcmm} estimation: model specification and goodness-of-fit statistics.
}
\usage{
\method{print}{hlme}(x,\dots)
\method{print}{Jointlcmm}(x,\dots)
\method{print}{lcmm}(x,\dots)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
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
\code{\link{hlme}}, \code{\link{Jointlcmm}}, \code{\link{lcmm}},
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{print}
