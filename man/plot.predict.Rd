\name{plot.predict}
\alias{plot.predict.hlme}
\alias{plot.predict.lcmm}
\alias{plot.predict.Jointlcmm}
\alias{plot.predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of class-specific marginal predictions for the longitudinal outcome 
}
\description{
This function provides the class-specific predicted trajectories stemmed from a \code{hlme}, \code{lcmm} or \code{Jointlcmm} object.
}
\usage{
\method{plot.predict}{hlme}(x,newdata,var.time,legend.loc="topright",ylim=NULL,\dots)
\method{plot.predict}{lcmm}(x,newdata,var.time,legend.loc="topright",ylim=NULL,\dots)
\method{plot.predict}{Jointlcmm}(x,newdata,var.time,legend.loc="topright",ylim=NULL,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object inheriting from classes \code{hlme}, \code{lcmm} or \code{Jointlcmm} representing respectively a fitted latent class linear mixed-effects model, a general latent class mixed model or a joint latent class mixed model.
}
  \item{newdata}{
data frame containing the data from which predictions are computed. Data frame should include all the covariates listed in x$Xnames (in the same order). Names of data frame should be exactly x$Xnames.
}
  \item{var.time}{
A character string containing the name of the variable that corresponds to time in the data frame (x axis in the plot).
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
\code{\link{hlme}}, \code{\link{lcmm}}, \code{\link{Jointlcmm}}   
}
\examples{

\dontrun{
################# Prediction from linear latent class model
data(data_hlme)
## fitted model
m<-lcmm(Y~Time+X1+X1_time,mixture=~Time,random=~Time,classmb=~X2+X3,
subject='ID',ng=2,data=data_hlme)
## newdata for predictions plot
newdata<-data.frame(intercept=rep(1,100),Time=seq(0,5,length=100),
X1=rep(0,100),X1_time=rep(0,100),X2=rep(0,100),X3=rep(0,100))
plot.predict(m,newdata,"Time","right")
## data from the first subject for predictions plot
firstdata<-cbind(intercept=rep(1,3),data_hlme[1:3,c("Time","X1",
"X1_time","X2","X3")])
plot.predict(m,firstdata,"Time","right")


################# Prediction from a joint latent class model
data(data_Jointlcmm)
## fitted model - see help of Jointlcmm function for details on the model
m3 <- Jointlcmm(fixed= Ydep1~Time+X1+X1_Time,mixture=~Time,random=~Time,
classmb=~X3,subject='ID',survival = Surv(Tevent,Event)~ mixture(X1)+X2
,hazard="5-equi-splines",hazardtype="PH",ng=3,data=data_Jointlcmm)
# class-specific predicted trajectories (with characteristics of subject ID=193)
data <- data_Jointlcmm[data_Jointlcmm$ID==193,]
plot.predict(m3,var.time="Time",newdata=data)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
