\name{plot.linkfunction}
\alias{plot.linkfunction.lcmm}
\alias{plot.linkfunction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of the estimated transformation between the outcome and the underlying latent process
}
\description{
Displays the estimated transformation specified in the option \code{link} of \code{lcmm} function. It corresponds to the (non)linear parameterized link assumed between the oberved longitudinal outcome and the underlying latent process and estimated as long as other regression parameters. 
}
\usage{
plot.linkfunction.lcmm(x,legend.loc="topright",...)
plot.linkfunction(x,legend.loc="topright",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object inheriting from class \code{lcmm}, representing a fitted latent class mixed model
}
  \item{legend.loc}{
keyword for the position of the legend from the list \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"},
     \code{"left"}, \code{"topleft"},\code{"top"}, \code{"topright"}, \code{"right"} and
     \code{"center"}. By default, the legend is located in the top right of the plot. 
}
  \item{\dots}{
other parameters to be passed through to plotting functions
}
}
\author{
Cecile Proust-Lima and Benoit Liquet
}
\seealso{
\code{\link{lcmm}}  
}

\examples{

\dontrun{


#### Estimation of homogeneous mixed models with different assumed link functions
#### quadratic mean trajectory and independent random intercept, slope and quadratic slope
#### (comparison of linear, Beta and 3 splines link functions)
# linear link function
m10<-lcmm(Ydep2~Time+Time_2,random=~Time+Time_2,subject='ID',ng=1,idiag=TRUE,data=data_Jointlcmm,link="linear")
summary(m10)
# Beta link function
m11<-lcmm(Ydep2~Time+Time_2,random=~Time+Time_2,subject='ID',ng=1,idiag=TRUE,data=data_Jointlcmm,link="beta")
summary(m11)
plot.linkfunction(m11)
# I-splines with 3 equidistant nodes
m12<-lcmm(Ydep2~Time+Time_2,random=~Time+Time_2,subject='ID',ng=1,idiag=TRUE,data=data_Jointlcmm,link="3-equi-splines")
summary(m12)
# I-splines with 5 nodes at quantiles
m13<-lcmm(Ydep2~Time+Time_2,random=~Time+Time_2,subject='ID',ng=1,idiag=TRUE,data=data_Jointlcmm,link="5-quant-splines")
summary(m13)
# I-splines with 5 nodes, and interior nodes entered manually
m14<-lcmm(Ydep2~Time+Time_2,random=~Time+Time_2,subject='ID',ng=1,idiag=TRUE,data=data_Jointlcmm,link="5-manual-splines",intnodes=c(10,20,25))
summary(m14)
plot.linkfunction(m14)


#### Plot of estimated different link functions:
#### (applicable for models that only differ in the "link function" used. 
#### Otherwise, the latent process scale is different and a rescaling is necessary)
transfo=data.frame(marker=m10$estimlink[,1],linear=m10$estimlink[,2],beta=m11$estimlink[,2]
,spl_3e=m12$estimlink[,2],spl_5q=m13$estimlink[,2],spl_5m=m14$estimlink[,2])
plot(transfo[,1]~transfo[,2],xlim=c(-10,5),col=1,type='l',xlab="latent process",ylab="marker")
par(new=TRUE)
plot(transfo[,1]~transfo[,3],xlim=c(-10,5),col=2,type='l',xlab="",ylab="")
par(new=TRUE)
plot(transfo[,1]~transfo[,4],xlim=c(-10,5),col=3,type='l',xlab="",ylab="")
par(new=TRUE)
plot(transfo[,1]~transfo[,5],xlim=c(-10,5),col=4,type='l',xlab="",ylab="")
legend(x="bottomright",legend=colnames(transfo[,2:5]),col=1:4,lty=1,inset=.02)





}
}
