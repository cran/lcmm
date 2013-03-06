\name{plot.baselinerisk}
\alias{plot.baselinerisk.Jointlcmm}
\alias{plot.baselinerisk}
\title{
Plot of the class-specific baseline risk functions estimated from a joint latent class mixed model
}
\description{
Displays the estimated baseline risk functions for the time-to-event of interest in each latent class
}
\usage{
\method{plot.baselinerisk}{Jointlcmm}(x,legend.loc="topleft",ylim=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object inheriting from class \code{Jointlcmm}, representing a fitted joint latent class mixed model
}
  \item{legend.loc}{
keyword for the position of the legend from the list \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"},
     \code{"left"}, \code{"topleft"},\code{"top"}, \code{"topright"}, \code{"right"} and
     \code{"center"}. By default, the legend is located in the top left of the plot. 
}
  \item{ylim}{
optional numeric vector of length 2, giving the y coordinate range.
}
  \item{\dots}{
other parameters to be passed through to plotting functions
}
}
\author{
Cecile Proust-Lima, Amadou Diakite and Benoit Liquet
}
\seealso{
\code{\link{Jointlcmm}}
}

\examples{

%\dontrun{
data(data_Jointlcmm)
#### estimation with 3 latent classes (ng=3) - see Jointlcmm 
#### help for details on the model
m3 <- Jointlcmm(fixed= Ydep1~Time*X1,mixture=~Time,random=~Time,
classmb=~X3,subject='ID',survival = Surv(Tevent,Event)~ X1+mixture(X2),
hazard="3-quant-splines",hazardtype="PH",ng=3,data=data_Jointlcmm,
B=c(0.7667, 0.4020, -0.8243, -0.2726, 0.0000, 0.0000, 0.0000, 0.3020, 
-0.6212, 2.6247, 5.3139, -0.0255, 1.3595, 0.8172, -11.6867, 10.1668,
10.2355, 11.5137, -2.6209, -0.4328, -0.6062, 1.4718, -0.0378, 0.8505,
0.0366, 0.2634, 1.4981))


# Class-specific predicted baseline risk & survival functions in the 
# 3-class model retained (for the reference value of the covariates) 
plot.baselinerisk(m3,bty="l")
plot.baselinerisk(m3,ylim=c(0,5),bty="l")
plot.survival(m3,bty="l")
%}
}




