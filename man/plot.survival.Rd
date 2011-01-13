\name{plot.survival}
\alias{plot.survival}
\alias{plot.survival.Jointlcmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of the class-specific event-free probabilities (survival functions) estimated from a joint latent class mixed model
}
\description{
Displays the estimated event-free probabilities (survival functions) for the time-to-event of interest in each latent class
}
\usage{
\method{plot.survival}{Jointlcmm}(x,legend.loc="topright",ylim=c(0,1),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{
an object inheriting from class \code{Jointlcmm}, representing a fitted joint latent class mixed model
}
  \item{legend.loc}{
keyword for the position of the legend from the list \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"},
     \code{"left"}, \code{"topleft"},\code{"top"}, \code{"topright"}, \code{"right"} and
     \code{"center"}. By default, the legend is located in the top right of the plot. 
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Jointlcmm}}
}

\examples{

\dontrun{

data(data_Jointlcmm)
#### estimation with 3 latent classes (ng=3) - see Jointlcmm help for details 
#### on the model
m3 <- Jointlcmm(fixed= Ydep1~Time+X1+X1_Time,mixture=~Time,random=~Time,
classmb=~X3,subject='ID',survival = Surv(Tevent,Event)~ mixture(X1)+X2,
hazard="5-equi-splines",hazardtype="PH",ng=3,data=data_Jointlcmm)
# Class-specific predicted baseline risk & survival functions in the 3-class 
# model retained (for the reference value of the covariates) 
plot.baselinerisk(m3)
plot.baselinerisk(m3,ylim=c(0,5))
plot.survival(m3)
}
}
