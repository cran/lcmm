\name{lcmm}
\alias{lcmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of latent class mixed-effect models for different types of outcomes (continuous Gaussian, continuous non-Gaussian or ordinal)
}
\description{
This function fits mixed models and latent class mixed models for different types of outcomes. It handles continuous longitudinal outcomes (Gaussian or non-Gaussian) as well as bounded quantitative, discrete and ordinal longitudinal outcomes. 
The different types of the outcomes are taken into account using parameterized nonlinear link functions between the observed outcome and the underlying latent process of interest it measures. 
At the latent process level, the model estimates a standard linear mixed model or a latent class mixed model when heterogeneity in the population is investigated (in the same way as in function \code{hlme}). 
Parameters of the nonlinear link function and of the latent process mixed model are estimated simultaneously using a maximum likelihood method. 
}
\usage{
lcmm(fixed, mixture, random, subject, classmb, ng = 1, idiag = FALSE, 
nwg = FALSE, link = "linear", intnodes = NULL, epsY = 0.5, data, B, 
convB = 1e-04, convL = 1e-04, convG = 1e-04, maxiter=500, nsim=100, prior)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fixed}{
a two-sided linear formula object for specifying the fixed-effects in the linear mixed model at the latent process level. The response outcome is on the left of \code{~} and the covariates are separated by \code{+} on the right of the \code{~}. 
Fo identifiability purposes, the intercept specified by default should not be removed by a \code{-1}.
}
  \item{mixture}{
a one-sided formula object for the class-specific fixed effects in the latent process mixed model (to specify only for a number of latent classes greater than 1). 
Among the list of covariates included in \code{fixed}, the covariates with class-specific regression parameters are entered in \code{mixture} separated by \code{+}. 
By default, an intercept is included. If no intercept, \code{-1} should be the first term included.
}
  \item{random}{
an optional one-sided formula for the random-effects in the latent process mixed model. Covariates with a random-effect are separated by \code{+}. 
By default, an intercept is included. If no intercept, \code{-1} should be the first term included.
}
  \item{subject}{name of the covariate representing the grouping structure.
}
  \item{classmb}{
an optional one-sided formula describing the covariates in the class-membership multinomial logistic model. Covariates included are separated by \code{+}. 
No intercept should be included in this formula.
}
  \item{ng}{
number of latent classes considered. If \code{ng=1} no \code{mixture} nor \code{classmb} should be specified. If \code{ng>1}, \code{mixture} is required.
}
  \item{idiag}{
optional logical for the variance-covariance structure of the random-effects. If \code{FALSE}, a non structured matrix of variance-covariance is considered (by default). 
If \code{TRUE} a diagonal matrix of variance-covariance is considered. 
}
 \item{nwg}{
optional logical of class-specific variance-covariance of the random-effects. If \code{FALSE} the variance-covariance matrix is common over latent classes (by default).
If \code{TRUE} a class-specific proportional parameter multiplies the variance-covariance matrix in each class (the proportional parameter in the last latent class equals 1 to ensure identifiability).
}
  \item{link}{
optional family of link functions to estimate. By default, "linear" option specifies a linear link function leading to a standard linear mixed model (homogeneous or heterogeneous as estimated in \code{hlme}).  Other possibilities include "beta" for estimating a link function from the family of Beta cumulative distribution functions, "thresholds" for using a threshold model to describe the correspondence between each level of an ordinal outcome and the underlying latent process, and "Splines" for approximating the link function by I-splines. For this latter case, the number of nodes and the nodes location should be also specified. The number of nodes is first entered followed by \code{-}, then the location is specified with "equi", "quant" or "manual" for respectively equidistant nodes, nodes at quantiles of the marker distribution or interior nodes entered manually in argument \code{intnodes}. It is followed by \code{-} and finally "splines" is indicated. For example, "7-equi-splines" means I-splines with 7 equidistant nodes, "6-quant-splines" means I-splines with 6 nodes located at the quantiles of the marker distribution and "9-manual-splines" means I-splines with 9 nodes, the vector of 7 interior nodes being entered in the argument \code{intnodes}.  
}
  \item{intnodes}{
optional vector of interior nodes. This argument is only required for a I-splines link function with nodes entered manually. 
}
  \item{epsY}{
optional definite positive real used to rescale the marker in (0,1) when the beta link function is used. By default, epsY=0.5.
}
  \item{data}{
optional data frame containing the variables named in \code{fixed}, \code{mixture}, \code{random}, \code{classmb} and \code{subject}. 
}
  \item{B}{
optional vector containing the initial values for the parameters. The order in which the parameters are included is detailed in \code{details} section. 
If no vector is specified, a preliminary analysis involving the estimation of a mixed model (\code{lcmm} with ng=1) is performed to choose initial values. 
Due to possible local maxima in latent class mixed models, the \code{B} vector should be specified and several different starting points should be tried when ng>1.
}

\item{convB}{optional threshold for the convergence criterion based on the parameter stability. By default, convB=0.0001
}
\item{convL}{optional threshold for the convergence criterion based on the log-likelihood stability. By default, convL=0.0001
}
\item{convG}{optional threshold for the convergence criterion based on the derivatives. By default, convG=0.0001
}
\item{maxiter}{optional maximum number of iterations for the Marquardt iterative algorithm. By default, maxiter=500
}
  \item{nsim}{
number of points used to plot the estimated link function. By default, nsim=100.
}
  \item{prior}{name of the covariate containing the prior on the latent class membership. The covariate should be an integer with values in 0,1,...,ng. When there is no prior, the value should be 0. When there is a prior for the subject, the value should be the number of the latent class (in 1,...,ng)
}

}

\details{

A. THE PARAMETERIZED LINK FUNCTIONS

\code{lcmm} function estimates latent class mixed models for different types of outcomes by assuming a parameterized link function for linking the outcome Y(t) with the underlying latent process L(t) it measures. To fix the latent process dimension, we chose to constrain the (first) intercept of the latent class mixed model at the latent process level at 0 and the standard error of the gaussian error of measurement at 1. These two parameters are replaced by additional parameters in the parameterized link function. 

1. With the "linear" link function, 2 parameters are required that correspond directly to the intercept and the standard error: (Y - b1)/b2 = L(t)

2. With the "beta" link function, 4 parameters are required for the following transformation: [ h(Y(t)',b1,b2) - b3]/b4 where h is the Beta CDF with canonical parameters c1 and c2 that can be derived from b1 and b2 as c1=exp(b1)/[exp(b2)*(1+exp(b1))] and c2=1/[exp(b2)*(1+exp(b1))], and Y(t)' is the rescaled outcome i.e. Y(t)'= [ Y(t) - min(Y(t)) + epsY ] / [ max(Y(t)) - min(Y(t)) +2*epsY ]   

3. With the "splines" link function, n+2 parameters are required for the following transformation b_1 + b_2*I_1(Y(t)) + ... + b_\{n+2\} I_\{n+1\}(Y(t)), where I_1,...,I_\{n+1\} is the basis of quadratic I-splines

4. With the "thresholds" link function for an ordinal outcome in levels 0,...,C, C parameters are required for the following transformation: Y(t)=c <=> b_{c} < L(t) <= b_\{c+1\} with b_0 = - infinity and b_\{C+1\}=+infinity. 

Details of these parameterized link functions can be found in the referred papers. 


B. THE VECTOR OF PARAMETERS B

The parameters in the vector of initial values \code{B} or in the vector of maximum likelihood estimates \code{best} are included in the following order: 
(1) ng-1 parameters are required for intercepts in the latent class membership model, and if covariates are included in \code{classmb}, ng-1 paramaters should be entered for each one; 
(2) for all covariates in \code{fixed}, one parameter is required if the covariate is not in \code{mixture}, ng paramaters are required if the covariate is also in \code{mixture}; When ng=1, the intercept is not estimated and no parameter should be specified in \code{B}. When ng>1, the first intercept is not estimated and only ng-1 parameters should be specified in \code{B};
(3) the variance of each random-effect specified in \code{random} (including the intercept) 
if \code{idiag=TRUE} and the inferior triangular variance-covariance matrix of all the random-effects if \code{idiag=FALSE}; 
(4) only if \code{nwg=TRUE}, ng-1 parameters for class-specific proportional coefficients
 for the variance covariance matrix of the random-effects; 
(5) In contrast with hlme, due to identifiability purposes, the standard error of the Gaussian error is not estimated (fixed at 1), and should not be specified in \code{B};
(6) The parameters of the link function: 2 for "linear", 4 for "beta", n+2 for "splines" with n nodes and the number of levels minus one for "thresholds".

We understand that it can be difficult to enter the correct number of parameters in \code{B} at the first place. So we recommend to run the program without specifying the initial vector \code{B} even if this model does not converge. As the final vector \code{best} has exactly the same structure as \code{B} (even when the program stops without convergence), it will help defining a satisfying vector of initial values \code{B} for next runs. 


C. CAUTIONS REGARDING THE USE OF THE PROGRAM

Some caution should be made when using the program. 
First, convergence criteria are very strict as they are based on derivatives of the log-likelihood in addition to the parameter and log-likelihood stability. 
In some cases, the program may not converge and reach the maximum number of iterations fixed at 100. 
In this case, the user should check that parameter estimates at the last iteration are not on the boundaries of the parameter space. 
If the parameters are on the boundaries of the parameter space, the identifiability of the model should be assessed. 
If not, the program should be run again with other initial values, with a higher maximum number of iterations or less strict convergence tolerances.

Specifically when investigating heterogeneity (that is with ng>1):
(1) As the log-likelihood of a latent class model can have multiple maxima, a careful choice of the initial values is crucial for ensuring convergence toward the global maximum. 
The program can be run without entering the vector of initial values (see point 2). 
However, we recommend to systematically enter initial values in \code{B} and try different sets of initial values. 
(2) The automatic choice of initial values we provide requires the estimation of a preliminary linear mixed model. The user should be aware that first, this preliminary analysis can take time for large datatsets and second, 
that the generated initial values can be very not likely and even may converge slowly to a local maximum. 
This is a reason why specification of initial values in \code{B} should be preferred.

}
\value{
The list returned is:
\item{ns}{number of grouping units in the dataset}
\item{ng}{number of latent classes}
\item{loglik}{log-likelihood of the model}
\item{best}{vector of parameter estimates in the same order as specified in \code{B} and detailed in section \code{details}}
\item{V}{vector containing the upper triangle matrix of variance-covariance estimates of \code{Best} with exception for variance-covariance parameters of the random-effects for which \code{V} contains the variance-covariance estimates of the Cholesky transformed parameters displayed in \code{cholesky}}
\item{gconv}{vector of convergence criteria: 1. on the parameters, 2. on the likelihood, 3. on the derivatives}
\item{conv}{status of convergence: =1 if the convergence criteria were satisfied, =2 if the maximum number of iterations was reached, =4 or 5 if a problem occured during optimisation}
\item{call}{the matched call}
\item{niter}{number of Marquardt iterations}
\item{dataset}{dataset}
\item{N}{internal information used in related functions}
\item{name.mat.cov}{internal information used in related functions}
\item{idiag}{internal information used in related functions}
\item{pred}{table of individual predictions and residuals; it includes marginal predictions (pred_m), marginal residuals (resid_m), subject-specific predictions (pred_ss) and subject-specific residuals 
(resid_ss) averaged over classes, the observation (obs) and finally the class-specific marginal and subject-specific predictions 
(with the number of the latent class: pred_m_1,pred_m_2,...,pred_ss_1,pred_ss_2,...)}
\item{pprob}{table of posterior classification and posterior individual class-membership probabilities}
\item{Xnames}{list of covariates included in the model - for use in function \code{\link{plot.predict.hlme}}}

\item{predRE}{table containing individual predictions of the random-effects : a column per random-effect, a line per subject}
\item{cholesky}{vector containing the estimates of the Cholesky transformed parameters of the variance-covariance matrix of the random-effects}
\item{estimlink}{table containing the simulated values of the marker and corresponding estimated link function}
\item{linktype}{indicator of link function type: 0 for linear, 1 for beta, 2 for splines and 3 for thresholds}
\item{linknodes}{vector of nodes useful only for the 'splines' link function}

}


\author{
Cecile Proust-Lima, Amadou Diakite and Benoit Liquet
}
%% ~Make other sections like Warning with \section{Warning }{....} ~


\references{
Proust C and Jacqmin-Gadda H. Estimation of linear mixed models with a mixture of distribution for the random-effects. Comput Methods Programs Biomed 78:165-73

Proust, C., Jacqmin-Gadda, H., Taylor, J. M., Ganiayre, J. and Commenges, D. (2006). A
nonlinear model with latent process for cognitive evolution using multivariate longitudinal
data. Biometrics 62, 1014-24.

Proust-Lima, C., Dartigues, J.-F. and Jacqmin-Gadda, H. (2011). Misuse of the linear mixed
model when evaluating risk factors of cognitive decline. Technical report
}


\seealso{

\code{\link{postprob}},\code{\link{plot.postprob}},\code{\link{plot.linkfunction}},\code{\link{plot.predict}},\code{\link{hlme}}
}
\examples{
\dontrun{
#### Estimation of homogeneous mixed models with different assumed link functions
#### quadratic mean trajectory and independent random intercept, slope and quadratic slope
#### (comparison of linear, Beta and 3 splines link functions)
data(data_Jointlcmm)
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


#### Estimation of 2-latent class mixed models with different assumed link functions
#### with individual and class specific linear trend
# Linear link function
m20<-lcmm(Ydep2~Time,random=~Time,subject='ID',mixture=~Time,ng=2,idiag=TRUE,data=data_Jointlcmm,link="linear")
summary(m20)
postprob(m20)
# Beta link function
m21<-lcmm(Ydep2~Time,random=~Time,subject='ID',mixture=~Time,ng=2,idiag=TRUE,data=data_Jointlcmm,link="beta")
summary(m21)
postprob(m21)
# I-splines link function (and 5 nodes at quantiles)
m22<-lcmm(Ydep2~Time,random=~Time,subject='ID',mixture=~Time,ng=2,idiag=TRUE,data=data_Jointlcmm,link="5-quant-splines")
summary(m22)
postprob(m22)

}
}
