\name{plot.linkfunction}
\alias{plot.linkfunction.lcmm}
\alias{plot.linkfunction}
\alias{plot.linkfunction.multlcmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of the estimated transformation between the outcome and the underlying latent process
}
\description{
Displays the estimated transformation specified in the option \code{link} of \code{lcmm} function. It corresponds to the (non)linear parameterized link assumed between the oberved longitudinal outcome and the underlying latent process and estimated as long as other regression parameters. 
}
\usage{
\method{plot.linkfunction}{lcmm}(x,legend.loc="topright",...)
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
Cecile Proust-Lima, Viviane Philipps and Benoit Liquet
}
\seealso{
\code{\link{lcmm}}  
}

\examples{

\dontrun{

#### Estimation of homogeneous mixed models with different assumed link
#### functions, a quadratic mean trajectory for the latent process with 
#### independent random intercept, slope and quadratic slope
#### (comparison of linear, Beta and 3 splines link functions)
data(data_Jointlcmm)
# linear link function
m10<-lcmm(Ydep2~Time+I(Time^2),random=~Time+I(Time^2),subject='ID',ng=1,
data=data_Jointlcmm,link="linear",
B=c(-0.7454, -0.2031,  0.2715,  0.2916 , 0.6114, -0.0064,  0.0545,
  0.0128, 25.3795, 2.2371))
summary(m10)
# Beta link function
m11<-lcmm(Ydep2~Time+I(Time^2),random=~Time+I(Time^2),subject='ID',ng=1,
data=data_Jointlcmm,link="beta",B=c(-0.9109, -0.0831,  0.5194,  0.1910 ,
 0.8984, -0.0179, -0.0636,  0.0045,  0.5514, -0.7692,  0.7037,  0.0899))
summary(m11)
plot.linkfunction(m11)
# I-splines with 3 equidistant nodes
m12<-lcmm(Ydep2~Time+I(Time^2),random=~Time+I(Time^2),subject='ID',ng=1,
data=data_Jointlcmm,link="3-equi-splines",B=c(-0.9272, -0.0753 , 0.5304, 
 0.1950,  0.9260, -0.0204, -0.0739 , 0.0059, -7.8369,  0.9228 ,-1.4689,
  2.0396,  1.8102))
summary(m12)
# I-splines with 5 nodes at quantiles
m13<-lcmm(Ydep2~Time+I(Time^2),random=~Time+I(Time^2),subject='ID',ng=1,
data=data_Jointlcmm,link="5-quant-splines",B=c(-0.9262, -0.0759,  0.5306,
  0.1927,  0.9364, -0.0191, -0.0764,  0.0063, -7.8144, -1.0691,  1.3714,
 1.9428,  1.2957,  0.9396,  1.0587))
summary(m13)
# I-splines with 5 nodes, and interior nodes entered manually
m14<-lcmm(Ydep2~Time+I(Time^2),random=~Time+I(Time^2),subject='ID',ng=1,
data=data_Jointlcmm,link="5-manual-splines",intnodes=c(10,20,25),
B=c(-0.9315, -0.0739 , 0.5254 , 0.1933,  0.9418, -0.0206, -0.0776,
  0.0064, -7.8645, 0.7470,  1.2080,  1.5537 , 1.7558 , 1.3386 , 1.0982))
summary(m14)
plot.linkfunction(m14,bty="l")



#### Plot of estimated different link functions:
#### (applicable for models that only differ in the "link function" used. 
#### Otherwise, the latent process scale is different and a rescaling
#### is necessary)
transfo <- data.frame(marker=m10$estimlink[,1],linear=m10$estimlink[,2],
beta=m11$estimlink[,2],spl_3e=m12$estimlink[,2],spl_5q=m13$estimlink[,2],
spl_5m=m14$estimlink[,2])
plot(transfo[,1]~transfo[,2],xlim=c(-10,5),col=1,type='l',
xlab="latent process",ylab="marker",bty="l")
par(new=TRUE)
plot(transfo[,1]~transfo[,3],xlim=c(-10,5),col=2,type='l',
xlab="",ylab="",bty="l")
par(new=TRUE)
plot(transfo[,1]~transfo[,4],xlim=c(-10,5),col=3,type='l',
xlab="",ylab="",bty="l")
par(new=TRUE)
plot(transfo[,1]~transfo[,5],xlim=c(-10,5),col=4,type='l',
xlab="",ylab="",bty="l")
legend(x="bottomright",legend=colnames(transfo[,2:5]),col=1:4,
lty=1,inset=.02,bty="n")

}
}
