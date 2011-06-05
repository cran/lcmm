\name{hlme}
\alias{hlme}
\title{
Estimation of latent class linear mixed models
}
\description{
This function fits a latent class linear mixed model (LCLMM) also known as growth mixture model or heterogeneous linear mixed model. 
The LCLMM consists in assuming that the population is divided in a finite number of latent classes; 
each latent class is characterised by a specific mean trajectory which described by a class-specific linear mixed model. 
Both the latent class membership and the trajectory can be explained according to covariates. 
This model is limited to a Gaussian outcome.
} 

\usage{
hlme(fixed, mixture, random, subject, classmb, ng = 1, idiag = FALSE, 
 nwg = FALSE, data, B,convB=0.0001,convL=0.0001,convG=0.0001,prior,maxiter=500)
}


\arguments{
  \item{fixed}{
a two-sided linear formula object for the fixed-effects in the linear mixed model. The response outcome is on the left of \code{~} and the covariates are separated by \code{+} on the right of the \code{~}. 
By default, an intercept is included. If no intercept, \code{-1} should be the first term included on the right of \code{~}.
}
  \item{mixture}{
a one-sided formula object for the class-specific fixed effects in the linear mixed model (to specify only for a number of latent classes greater than 1). 
Among the list of covariates included in \code{fixed}, the covariates with class-specific regression parameters are entered in \code{mixture} separated by \code{+}. 
By default, an intercept is included. If no intercept, \code{-1} should be the first term included.
}
  \item{random}{
an optional one-sided formula for the random-effects in the linear mixed model. Covariates with a random-effect are separated by \code{+}. 
By default, an intercept is included. If no intercept, \code{-1} should be the first term included.
}
  \item{subject}{name of the covariate representing the grouping structure
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
  \item{data}{
optional data frame containing the variables named in \code{fixed}, \code{mixture}, \code{random}, \code{classmb} and \code{subject}. 
}
  \item{B}{
an optional vector containing the initial values for the parameters. The order in which the parameters are included is detailed in \code{details} section. 
If no vector is specified, a preliminary analysis involving the estimation of a standard linear mixed model is performed to choose initial values. 
Due to possible local maxima, the \code{B} vector should be specified and several different starting points should be tried.
}
\item{convB}{optional threshold for the convergence criterion based on the parameter stability. By default, convB=0.0001
}
\item{convL}{optional threshold for the convergence criterion based on the log-likelihood stability. By default, convL=0.0001
}
\item{convG}{optional threshold for the convergence criterion based on the derivatives. By default, convG=0.0001
}
\item{prior}{name of the covariate containing the prior on the latent class membership. The covariate should be an integer with values in 0,1,...,ng. When there is no prior, the value should be 0. When there is a prior for the subject, the value should be the number of the latent class (in 1,...,ng)
}
\item{maxiter}{optional maximum number of iterations for the Marquardt iterative algorithm. By default, maxiter=300}
}




\details{
The parameters in the vector of initial values \code{B} or in the vector of maximum likelihood estimates \code{best} are included in the following order: 
(1) ng-1 parameters are required for intercepts in the latent class membership model, and if covariates are included in \code{classmb}, ng-1 paramaters should be entered for each one; 
(2) for all covariates in \code{fixed}, one parameter is required if the covariate is not in \code{mixture}, ng paramaters are required if the covariate is also in \code{mixture}; 
(3) the variance of each random-effect specified in \code{random} (including the intercept) 
if \code{idiag=TRUE} and the inferior triangular variance-covariance matrix of all the random-effects if \code{idiag=FALSE}; 
(4) only if \code{nwg=TRUE}, ng-1 parameters for class-specific proportional coefficients
 for the variance covariance matrix of the random-effects; 
(5) the variance of the residual error. 
We understand that it can be difficult to enter the correct number of parameters in \code{B} at the first place. So we recommend to run the program without specifying the initial vector \code{B} even if this model does not converge. As the final vector \code{best} has exactly the same structure as \code{B} (even when the program stops without convergence), it will help defining a satisfying vector of initial values \code{B} for next runs. 


Some caution should be made when using the program:


(1) As the log-likelihood of a latent class model can have multiple maxima, a careful choice of the initial values is crucial for ensuring convergence toward the global maximum. 
The program can be run without entering the vector of initial values (see point 2). 
However, we recommend to systematically enter initial values in \code{B} and try different sets of initial values. 


(2) The automatic choice of initial values we provide requires the estimation of a preliminary linear mixed model. The user should be aware that first, this preliminary analysis can take time for large datatsets and second, 
that the generated initial values can be very not likely and even may converge slowly to a local maximum. 
This is a reason why specification of initial values in \code{B} should be preferred.


(4) Convergence criteria are very strict as they are based on derivatives of the log-likelihood in addition to the parameter and log-likelihood stability. 
In some cases, the program may not converge and reach the maximum number of iterations fixed at 100. 
In this case, the user should check that parameter estimates at the last iteration are not on the boundaries of the parameter space. 
If the parameters are on the boundaries of the parameter space, the identifiability of the model should be assessed. 
If not, the program should be run again with other initial values, with a higher maximum number of iterations or less strict convergence tolerances.
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

}
\references{
Verbeke G and Lesaffre E (1996). A linear mixed-effects model with heterogeneity in
the random-effects population. Journal of the American Statistical Association 91, 217-21

Muthen B and Shedden K (1999). Finite mixture modeling with mixture outcomes using the EM algorithm. Biometrics 55, 463-9


Proust C and Jacqmin-Gadda H (2005). Estimation of linear mixed models with a mixture of
distribution for the random-effects. Computer Methods Programs Biomedicine 78, 165-73

}
\author{
Cecile Proust-Lima and Benoit Liquet
}
\seealso{
\code{\link{postprob}},\code{\link{plot.postprob}},\code{\link{summary}},\code{\link{plot.predict}}
}
\examples{
data(data_hlme)
### homogeneous linear mixed model (standard linear mixed model) 
### with independent random-effects
m1<-hlme(Y~Time+X1+X1_time,random=~Time,subject='ID',ng=1,idiag=TRUE,data=data_hlme)
m1
### latent class linear mixed model with 2 classes
m2<-hlme(Y~Time+X1+X1_time,mixture=~Time,random=~Time,classmb=~X2+X3,subject='ID',ng=2,data=data_hlme)
m2
### same model as m2 with the vector if initial values specified
m3<-hlme(Y~Time+X1+X1_time,mixture=~Time,random=~Time,classmb=~X2+X3,subject='ID',ng=2,data=data_hlme,B=c(0,0,0,30,25,0,-1,0,0,5,0,1,1))
m3
}


