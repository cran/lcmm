\name{plot.lcmm}
\alias{plot.hlme}
\alias{plot.lcmm}
\alias{plot.Jointlcmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots of marginal and subject-specific residuals 
}
\description{
This function displays 2 plots: the marginal residuals against the marginal predictions or the subject-specific residuals against the subject-specific predictions
}
\usage{
\method{plot}{hlme}(x,\dots)
\method{plot}{lcmm}(x,\dots)
\method{plot}{Jointlcmm}(x,\dots)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{x}{an object inheriting from classes \code{hlme}, \code{lcmm} or \code{Jointlcmm}, representing respectively a fitted latent class linear mixed model, a more general latent class mixed model or a joint latent class model
}
  \item{\dots}{
other parameters to be passed through to plotting functions
}
}
\value{
Returns two plots
}

\author{
Cecile Proust-Lima and Benoit Liquet
}



\seealso{
\code{\link{hlme}}, \code{\link{lcmm}}, \code{\link{Jointlcmm}}
}


\examples{
data(data_hlme)
m<-lcmm(Y~Time*X1,mixture=~Time,random=~Time,classmb=~X2+X3,
subject='ID',ng=2,data=data_hlme)
plot(m)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

