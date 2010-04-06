\name{plot.predict.hlme}
\alias{plot.predict.hlme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Marginal class-specific predictions plot
}
\description{
This function provides the class-specific predicted trajectories stemmed from a \code{hlme} object.
}
\usage{
\method{plot.predict}{hlme}(x,newdata,var.time,legend.loc="topright",\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object inheriting from class \code{hlme}, representing a fitted latent class
          linear mixed-effects model
}
  \item{newdata}{
data frame containing the data from which predictions are computed. Data frame should include all the covariates listed in x$Xnames (in the same order). Names of data frame should be exactly x$Xnames.
}
  \item{var.time}{
A character string containing the name of the variable that corresponds to time in the data frame (x axis in the plot)
}
  \item{legend.loc}{
keyword for the position of the legend from the list \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"},
     \code{"left"}, \code{"topleft"},\code{"top"}, \code{"topright"}, \code{"right"} and
     \code{"center"}. By default, the legend is located in the top right of the plot. 
}
  \item{\dots}{
further arguments to be passed to or from other methods.
          They are ignored in this function.
}
}
\value{
Returns a plot}
\author{
Cecile Proust-Lima and Benoit Liquet
}
\seealso{
\code{\link{hlme}}, \code{\link{postprob.hlme}}, \code{\link{plot.postprob.hlme}}   
}
\examples{
data(data_hlme)
## fitted model
m<-hlme(Y~Time+X1+X1_time,mixture=~Time,random=~Time,classmb=~X2+X3,
 subject=ID,ng=2,data=data_hlme,B=c(0,0,0,30,25,0,-1,0,0,5,0,1,1))
## newdata for predictions plot
newdata<-data.frame(intercept=rep(1,100),Time=seq(0,5,length=100),
 X1=rep(0,100),X1_time=rep(0,100),X2=rep(0,100),X3=rep(0,100))
plot.predict.hlme(m,newdata,"Time","right")
## data from the first subject for predictions plot
firstdata<-cbind(intercept=rep(1,3),
 data_hlme[1:3,c("Time","X1","X1_time","X2","X3")])
plot.predict.hlme(m,firstdata,"Time","right")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
