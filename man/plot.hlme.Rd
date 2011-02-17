\name{plot.hlme}
\alias{plot.hlme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots of marginal and subject-specific residuals 
}
\description{
This function displays 2 plots: the marginal residuals against the marginal predictions or the subject-specific residuals against the subject-specific predictions
}
\usage{
\method{plot}{hlme}(x,\dots)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object inheriting from class \code{hlme}, representing a fitted latent class
          linear mixed-effects model
}
  \item{\dots}{
other parameters to be passed through to plotting functions
}
}
\value{
Returns two plots}
\author{
Cecile Proust-Lima and Benoit Liquet
}

\seealso{
\code{\link{hlme}}
}
\examples{
data(data_hlme)
m<-hlme(Y~Time+X1+X1_time,mixture=~Time,random=~Time,classmb=~X2+X3,subject=ID,ng=2,data=data_hlme,B=c(0,0,0,30,25,0,-1,0,0,5,0,1,1))

plot(m)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

