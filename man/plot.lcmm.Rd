\name{plot.lcmm}
\alias{plot.hlme}
\alias{plot.lcmm}
\alias{plot.Jointlcmm}
\alias{plot.multlcmm}
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
\method{plot}{multlcmm}(x,\dots)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{x}{an object inheriting from classes \code{hlme}, \code{lcmm}, \code{multlcmm} or \code{Jointlcmm}, representing respectively a fitted latent class linear mixed model, a more general latent class mixed model or a joint latent class model
}
  \item{\dots}{
other parameters to be passed through to plotting functions
}
}
\value{
Returns two plots
}

\author{
Cecile Proust-Lima, Viviane Philipps and Benoit Liquet
}



\seealso{
\code{\link{hlme}}, \code{\link{lcmm}}, \code{\link{multlcmm}}, \code{\link{Jointlcmm}}
}


\examples{
data(data_hlme)
m<-lcmm(Y~Time*X1,mixture=~Time,random=~Time,classmb=~X2+X3,
subject='ID',ng=2,data=data_hlme,B=c(0.41,0.55,-0.18,-0.41,
-14.26,-0.34,1.33,13.51,24.65,2.98,1.18,26.26,0.97))
plot(m)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

