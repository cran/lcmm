\name{summary.lcmm}
\alias{summary.hlme}
\alias{summary.lcmm}
\alias{summary.Jointlcmm}
\alias{summary.epoce}
\alias{summary.Diffepoce}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of a \code{hlme}, \code{lcmm}, \code{Jointlcmm} \code{epoce} or \code{Diffepoce} objects
}
\description{
The function provides a summary of \code{hlme}, \code{lcmm} and \code{Jointlcmm} estimations, or \code{epoce} and \code{Diffepoce} computations.
}
\usage{
\method{summary}{hlme}(object,\dots)
\method{summary}{lcmm}(object,\dots)
\method{summary}{Jointlcmm}(object,\dots)
\method{summary}{epoce}(object,\dots)
\method{summary}{Diffepoce}(object,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object inheriting from classes \code{hlme}, \code{lcmm} for fitted latent class mixed-effects, or class \code{Jointlcmm}for a Joint latent class mixed model or \code{epoce} or \code{Diffepoce} for predictive accuracy computations.
}
  \item{\dots}{
further arguments to be passed to or from other methods.
          They are ignored in this function.
}
}
\author{
Cecile Proust-Lima, Amadou Diakite and Benoit Liquet
}
\seealso{
\code{\link{hlme}}, \code{\link{lcmm}}, \code{\link{Jointlcmm}}, \code{epoce}, \code{Diffepoce}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{print}
