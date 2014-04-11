\name{hlme}
\alias{hlme}
\title{
Estimation of latent class linear mixed models
}
\description{
This function fits linear mixed models and latent class linear mixed models (LCLMM) also known as growth mixture models or heterogeneous linear mixed models. 
The LCLMM consists in assuming that the population is divided in a finite number of latent classes. Each latent class is characterised by a specific trajectory modelled by a class-specific linear mixed model. 
Both the latent class membership and the trajectory can be explained according to covariates. 
This function is limited to a mixture of Gaussian outcomes. For other types of outcomes, please see function \code{lcmm}. For multivariate longitudinal outcomes, please see \code{multlcmm}.
} 

\usage{
hlme(fixed, mixture, random, subject, classmb, ng = 1,
 idiag = FALSE, nwg = FALSE, cor=NULL, data, B,
 convB=0.0001,convL=0.0001,convG=0.0001,prior,
 maxiter=500, subset=NULL, na.action=1)
}


\arguments{
  \item{fixed}{
two-sided linear formula object for the fixed-effects in the linear mixed model. The response outcome is on the left of \code{~} and the covariates are separated by \code{+} on the right of \code{~}. 
By default, an intercept is included. If no intercept, \code{-1} should be the first term included on the right of \code{~}.
}
  \item{mixture}{
one-sided formula object for the class-specific fixed effects in the linear mixed model (to specify only for a number of latent classes greater than 1). 
Among the list of covariates included in \code{fixed}, the covariates with class-specific regression parameters are entered in \code{mixture} separated by \code{+}. 
By default, an intercept is included. If no intercept, \code{-1} should be the first term included.
}
  \item{random}{
optional one-sided formula for the random-effects in the linear mixed model. Covariates with a random-effect are separated by \code{+}. 
By default, an intercept is included. If no intercept, \code{-1} should be the first term included.
}
  \item{subject}{name of the covariate representing the grouping structure specified with ''.
}
  \item{classmb}{
optional one-sided formula describing the covariates in the class-membership multinomial logistic model. Covariates included are separated by \code{+}. 
No intercept should be included in this formula.
}
  \item{ng}{
optional number of latent classes considered. If \code{ng=1} (by default) no \code{mixture} nor \code{classmb} should be specified. If \code{ng>1}, \code{mixture} is required.
}
  \item{idiag}{
optional logical for the structure of the variance-covariance matrix of the random-effects. If \code{FALSE}, a non structured matrix of variance-covariance is considered (by default). 
If \code{TRUE} a diagonal matrix of variance-covariance is considered. 
}
 \item{nwg}{
optional logical indicating if the variance-covariance of the random-effects is class-specific. If \code{FALSE} the variance-covariance matrix is common over latent classes (by default).
If \code{TRUE} a class-specific proportional parameter multiplies the variance-covariance matrix in each class (the proportional parameter in the last latent class equals 1 to ensure identifiability).
}
  \item{cor}{
  optional brownian motion or autoregressive process modeling the correlation between the observations.
  "BM" or "AR" should be specified, followed by the time variable between brackets. By default, no correlation is added. 
}
  \item{data}{
optional data frame containing the variables named in \code{fixed}, \code{mixture}, \code{random}, \code{classmb} and \code{subject}. 
}
  \item{B}{
optional vector containing the initial values for the parameters. The order in which the parameters are included is detailed in \code{details} section. 
If no vector is specified, a preliminary analysis involving the estimation of a standard linear mixed model is performed to choose initial values. 
Due to possible local maxima, the \code{B} vector should be specified and several different starting points should be tried.
}
\item{convB}{optional threshold for the convergence criterion based on the parameter stability. By default, convB=0.0001.
}
\item{convL}{optional threshold for the convergence criterion based on the log-likelihood stability. By default, convL=0.0001.
}
\item{convG}{optional threshold for the convergence criterion based on the derivatives. By default, convG=0.0001.
}
\item{prior}{optional name of a covariate containing a prior information about the latent class membership. The covariate should be an integer with values in 0,1,...,ng. Value 0 indicates no prior for the subject while a value in 1,...,ng indicates that the subject belongs to the corresponding latent class. 
}
\item{maxiter}{optional maximum number of iterations for the Marquardt iterative algorithm. By default, maxiter=500.}
  \item{subset}{
a specification of the rows to be used: defaults to all rows. This can be any valid indexing vector for the rows of data or if that is not supplied, a data frame made up of the variable used in formula.}
  \item{na.action}{
Integer indicating how NAs are managed. The default is 1 for 'na.omit'. The alternative is 2 for 'na.fail'. Other options such as 'na.pass' or 'na.exclude' are not implemented in the current version.}
}




\details{

A. THE VECTOR OF PARAMETERS B

The parameters in the vector of initial values \code{B} or equivalently in the vector of maximum likelihood estimates \code{best} are included in the following order: 

(1) ng-1 parameters are required for intercepts in the latent class membership model, and when covariates are included in \code{classmb}, ng-1 paramaters should be entered for each covariate; 

(2) for all covariates in \code{fixed}, one parameter is required if the covariate is not in \code{mixture}, ng paramaters are required if the covariate is also in \code{mixture}; 

(3) the variance of each random-effect specified in \code{random} (including the intercept) 
when \code{idiag=TRUE}, or the inferior triangular variance-covariance matrix of all the random-effects when \code{idiag=FALSE}; 

(4) only when \code{nwg=TRUE}, ng-1 parameters are required for the ng-1 class-specific proportional coefficients
 in the variance covariance matrix of the random-effects; 

(5) when \code{cor} is specified, 1 parameter corresponding to the variance of the Brownian motion should be entered with \code{cor=BM} and 2 parameters corresponding to the correlation and the variance parameters of the autoregressive process should be entered 

(6) the standard error of the residual error. 

We understand that it can be difficult to enter the correct number of parameters in \code{B} at the first place. So we recommend to run the program without specifying the initial vector \code{B} even if this model does not converge. To save time, the number of iterations can be restricted to 0 with \code{maxiter=0}. As the final vector \code{best} has exactly the same structure as \code{B} (even when the program stops without convergence), it will help defining a satisfying vector of initial values \code{B} for next runs. 


B. CAUTIONS

Some caution should be made when using the program:


(1) As the log-likelihood of a latent class model can have multiple maxima, a careful choice of the initial values is crucial for ensuring convergence toward the global maximum. 
The program can be run without entering the vector of initial values (see point 2). 
However, we recommend to systematically enter initial values in \code{B} and try different sets of initial values. 


(2) The automatic choice of initial values we provide requires the estimation of a preliminary linear mixed model. The user should be aware that first, this preliminary analysis can take time for large datatsets and second, 
that the generated initial values can be very not likely and even may converge slowly to a local maximum. 
This is a reason why the specification of initial values in \code{B} should be systematically preferred.


(4) Convergence criteria are very strict as they are based on the derivatives of the log-likelihood in addition to the parameter stability and log-likelihood stability. 
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
\item{idiag}{internal information used in related functions}
\item{pred}{table of individual predictions and residuals; it includes marginal predictions (pred_m), marginal residuals (resid_m), subject-specific predictions (pred_ss) and subject-specific residuals 
(resid_ss) averaged over classes, the observation (obs) and finally the class-specific marginal and subject-specific predictions 
(with the number of the latent class: pred_m_1,pred_m_2,...,pred_ss_1,pred_ss_2,...)}
\item{pprob}{table of posterior classification and posterior individual class-membership probabilities}
\item{Xnames}{list of covariates included in the model}

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
Cecile Proust-Lima, Benoit Liquet and Viviane Philipps
}
\seealso{
\code{\link{postprob}}, \code{\link{plot.postprob}}, \code{\link{summary}}, \code{\link{plot.predict}}
}
\examples{
%\dontrun{
##### Example of a latent class model estimated for a varying number
# of latent classes: 
# The model includes a subject- (ID) and class-specific linear 
# trend (intercept and Time in fixed, random and mixture components)
# and a common effect of X1 and its interaction with time over classes 
# (in fixed). 
# The random-effects are assumed independent (idiag=TRUE). The variance
# of the random intercept and slope are assumed to be equal 
# over classes (nwg=F).
# The covariate X3 predicts the class membership (in classmb). 
# !CAUTION: for illustration, only default initial values where used but 
# other sets of initial values should be tried to ensure convergence
# towards the global maximum.

data(data_hlme)

### homogeneous linear mixed model (standard linear mixed model) 
### with independent random-effects
m1<-hlme(Y~Time*X1,random=~Time,subject='ID',ng=1,idiag=TRUE,
data=data_hlme)
summary(m1)
### latent class linear mixed model with 2 classes
m2<-hlme(Y~Time*X1,mixture=~Time,random=~Time,classmb=~X2+X3,
subject='ID',ng=2,data=data_hlme,B=c(0.11,-0.74,-0.07,20.71,
29.39,-1,0.13,2.45,-0.29,4.5,0.36,0.79,0.97))
m2
summary(m2)
postprob(m2)
### same model as m2 with the vector of initial values specified
m3<-hlme(Y~Time*X1,mixture=~Time,random=~Time,classmb=~X2+X3,
subject='ID',ng=2,data=data_hlme,B=c(0,0,0,30,25,0,-1,0,0,5,0,1,1))
m3
summary(m3)

## plot of predicted trajectories using some newdata
newdata<-data.frame(Time=seq(0,5,length=100),
X1=rep(0,100),X2=rep(0,100),X3=rep(0,100))
plot.predict(m3,newdata,"Time","right",bty="l")

%}
}
