\name{plot.predict}
\alias{plot.predict.hlme}
\alias{plot.predict.lcmm}
\alias{plot.predict.Jointlcmm}
\alias{plot.predict.multlcmm}
\alias{plot.predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of class-specific marginal predictions for the longitudinal outcome 
}
\description{
This function provides the class-specific predicted trajectories stemmed from a \code{hlme}, \code{lcmm}, \code{multlcmm} or \code{Jointlcmm} object.
}
\usage{
\method{plot.predict}{hlme}(x,newdata,var.time,legend.loc="topright",
na.action=1,legend,add=FALSE,confint=FALSE,...)
\method{plot.predict}{lcmm}(x,newdata,var.time,legend.loc="topright",
na.action=1,legend,add=FALSE,confint=FALSE,...)
\method{plot.predict}{Jointlcmm}(x,newdata,var.time,legend.loc="topright",
na.action=1,legend,add=FALSE,confint=FALSE,...)
\method{plot.predict}{multlcmm}(x,newdata,var.time,legend.loc="topright",
na.action=1,legend,add=FALSE,confint=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object inheriting from classes \code{hlme}, \code{lcmm} or \code{Jointlcmm} representing respectively a fitted latent class linear mixed-effects model, a general latent class mixed model or a joint latent class mixed model.
}
  \item{newdata}{
data frame containing the data from which predictions are computed. The data frame should include at least all the covariates listed in x$Xnames2. Names in the data frame should be exactly x$Xnames2 that are the names of covariates specified in \code{lcmm}, \code{hlme}, \code{Jointlcmm} or \code{multlcmm} calls.
}
  \item{var.time}{
A character string containing the name of the variable that corresponds to time in the data frame (x axis in the plot).
}
  \item{legend.loc}{
keyword for the position of the legend from the list \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"},
     \code{"left"}, \code{"topleft"},\code{"top"}, \code{"topright"}, \code{"right"} and
     \code{"center"}. By default, the legend is located in the top right of the plot. 
}
 \item{na.action}{
Integer indicating how NAs are managed. The default is 1 for 'na.omit'. The alternative is 2 for 'na.fail'. Other options such as 'na.pass' or 'na.exclude' are not implemented in the current version.
}
  \item{legend}{
character or expression to appear in the legend. If no legend should be added, \code{"legend"} should be NULL.
}
  \item{add}{
logical indicating if the curves should be added to an existing plot. Default to FALSE.
}
  \item{confint}{
logical indicating if confidence intervals for the predictions should be calculated and plotted. Default to FALSE.
}
  \item{\dots}{
other parameters to be passed through to plotting functions or to legend 
}
}
\value{
Returns a data frame containing the time values at which the predictions are computated (corresponding to the x-axis in the plot), 
the predicted values of the latent process, the lower and the upper limits of the confidence intervals (if calculated).
}
\author{
Cecile Proust-Lima, Benoit Liquet and Viviane Philipps 
}
\seealso{
\code{\link{hlme}}, \code{\link{lcmm}}, \code{\link{Jointlcmm}}, \code{\link{multlcmm}}   
}
\examples{
% \dontrun{
################# Prediction from linear latent class model
data(data_hlme)
## fitted model
m<-lcmm(Y~Time*X1,mixture=~Time,random=~Time,classmb=~X2+X3,
subject='ID',ng=2,data=data_hlme,B=c(0.41,0.55,-0.18,-0.41,
-14.26,-0.34,1.33,13.51,24.65,2.98,1.18,26.26,0.97))
## newdata for predictions plot
newdata<-data.frame(Time=seq(0,5,length=100),
X1=rep(0,100),X2=rep(0,100),X3=rep(0,100))
plot.predict(m,newdata,"Time","right",bty="l")
## data from the first subject for predictions plot
firstdata<-data_hlme[1:3,]
plot.predict(m,firstdata,"Time","right",bty="l")


################# Prediction from a joint latent class model
data(data_Jointlcmm)
## fitted model - see help of Jointlcmm function for details on the model
m3 <- Jointlcmm(fixed= Ydep1~Time*X1,mixture=~Time,random=~Time,
classmb=~X3,subject='ID',survival = Surv(Tevent,Event)~X1+mixture(X2),
hazard="3-quant-splines",hazardtype="PH",ng=3,data=data_Jointlcmm,
B=c(0.7667, 0.4020, -0.8243, -0.2726, 0.0000, 0.0000, 0.0000, 0.3020,
-0.6212, 2.6247, 5.3139, -0.0255, 1.3595, 0.8172, -11.6867, 10.1668,
10.2355, 11.5137, -2.6209, -0.4328, -0.6062, 1.4718, -0.0378, 0.8505,
0.0366, 0.2634, 1.4981))
# class-specific predicted trajectories 
#(with characteristics of subject ID=193)
data <- data_Jointlcmm[data_Jointlcmm$ID==193,]
plot.predict(m3,var.time="Time",newdata=data,bty="l")
%}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
