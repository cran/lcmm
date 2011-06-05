\name{postprob}
\alias{postprob}
\alias{postprob.hlme}
\alias{postprob.lcmm}
\alias{postprob.Jointlcmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior classification stemmed from a \code{hlme}, \code{lcmm} or \code{Jointlcmm} estimation
}
\description{
This function provides informations about the posterior classification stemmed from a \code{hlme}, \code{lcmm} or \code{Jointlcmm} object.
}
\usage{
postprob(x,\dots)
postprob.hlme(x,\dots)
postprob.lcmm(x,\dots)
postprob.Jointlcmm(x,\dots)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{x}{an object inheriting from classes \code{hlme}, \code{lcmm} or \code{Jointlcmm} representing respectively a fitted latent class
          linear mixed-effects model, a more general latent class mixed model or a joint latent class model.}
 \item{\dots}{
further arguments to be passed to or from other methods.
          They are ignored in this function.
}
}

\details{
This function provides the number of subjects classified a posteriori in each latent class, and the classification table that contains the mean of the posterior probability of belonging to each latent class over the subjects classified in each of the latent classes. This table aims at evaluating the quality of the posterior classification. 
For \code{hlme}, \code{lcmm} objects, the posterior classification and the classification table are derived from the posterior class-membership probabilities given the vector of repeated measures that are contained in pprob output matrix.
For a \code{Jointlcmm} object, the first posterior classification and the classification table are derived from the posterior class-membership probabilities given the vector of repeated measures and the time-to-event information (that are contained in columns probYT1, probYT2, etc in pprob output matrix). The second posterior classification is derived from the posterior class-membership probabilities given only the vector of repeated measures (that are contained in columns probY1, probY2, etc in pprob output matrix).


}
\author{
Cecile Proust-Lima and Benoit Liquet
}
\note{
This function can only be used with latent class mixed models and joint latent class mixed models that include at least 2 latent classes
}
\seealso{
\code{\link{Jointlcmm}}, \code{\link{lcmm}}, \code{\link{hlme}},\code{\link{plot.postprob}}  
}
\examples{
data(data_hlme)
m<-lcmm(Y~Time+X1+X1_time,mixture=~Time,random=~Time,classmb=~X2+X3,subject='ID',ng=2,data=data_hlme)
postprob(m)
}
