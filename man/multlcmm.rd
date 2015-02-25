\name{multlcmm}
\alias{multlcmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of mutlivariate mixed-effect models and multivariate latent class mixed-effect models for multivariate longitudinal outcomes of possibly multiple types (continuous Gaussian, continuous non-Gaussian - curvilinear) that measure the same underlying latent process.
}
\description{
This function constitutes a multivariate extension of function \code{lcmm}. It fits multivariate mixed models and multivariate latent class mixed models for multivariate longitudinal outcomes of different types. It handles continuous longitudinal outcomes (Gaussian or non-Gaussian, curvilinear) as well as bounded quantitative and discrete longitudinal outcomes. Next version will also handle ordinal outcomes. 
The model assumes that all the outcomes measure the same underlying latent process defined as their common factor, and each outcome is related to this latent common factor by a specific parameterized link function. 
At the latent process level, the model estimates a standard linear mixed model or a latent class linear mixed model when heterogeneity in the population is investigated (in the same way as in function \code{hlme}).
Parameters of the nonlinear link functions and of the latent process mixed model are estimated simultaneously using a maximum likelihood method.
}
\usage{
multlcmm(fixed, mixture, random, subject, classmb, ng = 1, 
idiag = FALSE,nwg = FALSE, randomY=FALSE, link = "linear", 
intnodes = NULL, epsY = 0.5, cor=NULL, data, B, convB = 1e-04, 
convL = 1e-04, convG = 1e-04, maxiter=100,
nsim=100, prior,range=NULL, subset=NULL, na.action=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fixed}{
a two-sided linear formula object for specifying the fixed-effects in the linear mixed model at the latent process level. The response outcomes are separated by \code{+} on the left of \code{~} and the covariates are separated by \code{+} on the right of the \code{~}.
For identifiability purposes, the intercept specified by default should not be removed by a \code{-1}. Variables on which a contrast above the different outcomes should also be estimated are included with \code{contrast()}.
}
  \item{mixture}{
a one-sided formula object for the class-specific fixed effects in the latent process mixed model (to specify only for a number of latent classes greater than 1).
Among the list of covariates included in \code{fixed}, the covariates with class-specific regression parameters are entered in \code{mixture} separated by \code{+}.
By default, an intercept is included. If no intercept, \code{-1} should be the first term included.
}
  \item{random}{
an optional one-sided formula for the random-effects in the latent process mixed model. At least one random effect should be included for identifiability purposes. Covariates with a random-effect are separated by \code{+}.
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
 \item{randomY}{
 optional logical for including an outcome-specific random intercept. If \code{FALSE} no outcome-specific random intercept is added (default). If \code{TRUE} independent outcome-specific random intercepts with parameterized variance are included.
 }
  \item{link}{
optional vector of families of parameterized link functions to estimate (one by outcome). Option "linear" (by default) specifies a linear link function. Other possibilities include "beta" for estimating a link function from the family of Beta cumulative distribution functions and "Splines" for approximating the link function by I-splines. For this latter case, the number of nodes and the nodes location should be also specified. The number of nodes is first entered followed by \code{-}, then the location is specified with "equi", "quant" or "manual" for respectively equidistant nodes, nodes at quantiles of the marker distribution or interior nodes entered manually in argument \code{intnodes}. It is followed by \code{-} and finally "splines" is indicated.
 For example, "7-equi-splines" means I-splines with 7 equidistant nodes, "6-quant-splines" means I-splines with 6 nodes located at the quantiles of the marker distribution and "9-manual-splines" means I-splines with 9 nodes, the vector of 7 interior nodes being entered in the argument \code{intnodes}.
}
  \item{intnodes}{
optional vector of interior nodes. This argument is only required for a I-splines link function with nodes entered manually.
}
  \item{epsY}{
optional definite positive real used to rescale the marker in (0,1) when the beta link function is used. By default, epsY=0.5.
}
  \item{cor}{
  optional indicator for inclusion of an autocorrelated Gaussian process in the latent process linear (latent process) mixed model. Option "BM" indicates a brownian motion with parameterized variance. Option "AR" specifies an autoregressive process of order 1 with parameterized variance and correlation intensity. Each option should be followed by the time variable in brackets as \code{cor=BM(time)}. By default, no autocorrelated Gaussian process is added.
}
  \item{data}{
 data frame containing the variables named in \code{fixed}, \code{mixture}, \code{random}, \code{classmb} and \code{subject}.
}
  \item{B}{
optional vector containing the initial values for the parameters. The order in which the parameters are included is detailed in \code{details} section.
If no vector is specified and ng>1, a preliminary analysis involving the estimation of a linear mixed model (\code{multlcmm} with ng=1) is performed to choose initial values.
Due to possible local maxima in latent class mixed models, the \code{B} vector should be specified and several different starting points should be tried when ng>1.
}

  \item{convB}{optional threshold for the convergence criterion based on the parameter stability. By default, convB=0.0001.
}
  \item{convL}{optional threshold for the convergence criterion based on the log-likelihood stability. By default, convL=0.0001.
}
  \item{convG}{optional threshold for the convergence criterion based on the derivatives. By default, convG=0.0001.
}
  \item{maxiter}{optional maximum number of iterations for the Marquardt iterative algorithm. By default, maxiter=100.
}
  \item{nsim}{
number of points used to plot the estimated link functions. By default, nsim=100.
}
  \item{prior}{name of the covariate containing the prior on the latent class membership. The covariate should be an integer with values in 0,1,...,ng. When there is no prior, the value should be 0. When there is a prior for the subject, the value should be the number of the latent class (in 1,...,ng).
}
  \item{range}{optional vector indicating the range of the outcomes (that is the minimum and maximum). By default, the range is defined according to the minimum and maximum observed values of the outcome. The option should be used only for Beta and Splines transformations.
  }
  \item{subset}{optional vector giving the subset of observations in \code{data} to use. By default, all lines.
}
  \item{na.action}{
Integer indicating how NAs are managed. The default is 1 for 'na.omit'. The alternative is 2 for 'na.fail'. Other options such as 'na.pass' or 'na.exclude' are not implemented in the current version.
}
}

\details{

A. THE PARAMETERIZED LINK FUNCTIONS

\code{multlcmm} function estimates multivariate latent class mixed models for different types of outcomes by assuming a parameterized link function for linking each outcome Y_k(t) with the underlying latent common factor L(t) they measure. To fix the latent process dimension, we chose to constrain  at the latent process level the (first) intercept of the latent class mixed model at 0 and the standard error of the first random effect  at 1. 

1. With the "linear" link function, 2 parameters are required for the following transformation (Y(t) - b1)/b2 

2. With the "beta" link function, 4 parameters are required for the following transformation: [ h(Y(t)',b1,b2) - b3]/b4 where h is the Beta CDF with canonical parameters c1 and c2 that can be derived from b1 and b2 as c1=exp(b1)/[exp(b2)*(1+exp(b1))] and c2=1/[exp(b2)*(1+exp(b1))], and Y(t)' is the rescaled outcome i.e. Y(t)'= [ Y(t) - min(Y(t)) + epsY ] / [ max(Y(t)) - min(Y(t)) +2*epsY ].

3. With the "splines" link function, n+2 parameters are required for the following transformation b_1 + b_2*I_1(Y(t)) + ... + b_\{n+2\} I_\{n+1\}(Y(t)), where I_1,...,I_\{n+1\} is the basis of quadratic I-splines.
To constraint the parameters to be positive, except for b_1, the program estimates b_k^* (for k=2,...,n+2) so that b_k=(b_k^*)^2. This parameterization may lead in some cases to problems of convergence that we are currently addressing.

Details of these parameterized link functions can be found in the papers: Proust-Lima et al. (Biometrics 2006) and Proust-Lima et al. (BJMSP 2013).


B. THE VECTOR OF PARAMETERS B

The parameters in the vector of initial values \code{B} or in the vector of maximum likelihood estimates \code{best} are included in the following order:
(1) ng-1 parameters are required for intercepts in the latent class membership model, and if covariates are included in \code{classmb}, ng-1 paramaters should be entered for each one;
(2) for all covariates in \code{fixed}, one parameter is required if the covariate is not in \code{mixture}, ng paramaters are required if the covariate is also in \code{mixture}; When ng=1, the intercept is not estimated and no parameter should be specified in \code{B}. When ng>1, the first intercept is not estimated and only ng-1 parameters should be specified in \code{B};
(3) for all covariates included with \code{contrast()} in \code{fixed}, one supplementary parameter per outcome is required excepted for the last outcome for which the parameter is not estimated but deduced from the others;
(4) if \code{idiag=TRUE}, the variance of each random-effect specified in \code{random} is required excepted the first one (usually the intercept) which is constrained to 1.
(5) if \code{idiag=FALSE}, the inferior triangular variance-covariance matrix of all the random-effects is required excepted the first variance (usually the intercept) which is constrained to 1.
(5) only if \code{nwg=TRUE} and \code{ng}>1, ng-1 parameters for class-specific proportional coefficients
 for the variance covariance matrix of the random-effects;
(6) if \code{cor} is specified, the standard error of the Brownian motion or the standard error and the correlation parameter of the autoregressive process;
(7) the standard error of the outcome-specific Gaussian errors (one per outcome);
(8) if \code{randomY=TRUE}, the standard error of the outcome-specific random intercept (one per outcome);
(9) the parameters of each parameterized link function: 2 for "linear", 4 for "beta", n+2 for "splines" with n nodes.

We understand that it can be difficult to enter the correct number of parameters in \code{B} at the first place. So we recommend to run the program without specifying the initial vector \code{B} even if this model does not converge (maxiter can be even at 1 to run only 1 iteration of the optimization program). As the final vector \code{best} has exactly the same structure as \code{B} (even when the program stops without convergence), it will help defining a satisfying vector of initial values \code{B} for next runs.


C. CAUTIONS REGARDING THE USE OF THE PROGRAM

Some caution should be made when using the program.
First, convergence criteria are very strict as they are based on the derivatives of the log-likelihood in addition to the parameter and log-likelihood stability.
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
This is a reason why specification of initial values in \code{B} should be favoured.

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
\item{N}{internal information used in related functions}
\item{idiag}{internal information used in related functions}
\item{pred}{table of individual predictions and residuals in the underlying latent process scale; it includes marginal predictions (pred_m), marginal residuals (resid_m), subject-specific predictions (pred_ss) and subject-specific residuals
(resid_ss) averaged over classes, the transformed observations in the latent process scale (obs) and finally the class-specific marginal and subject-specific predictions
(with the number of the latent class: pred_m_1,pred_m_2,...,pred_ss_1,pred_ss_2,...).}
\item{pprob}{table of posterior classification and posterior individual class-membership probabilities}
\item{Xnames}{list of covariates included in the model}
\item{predRE}{table containing individual predictions of the random-effects : a column per random-effect, a line per subject.}
\item{cholesky}{vector containing the estimates of the Cholesky transformed parameters of the variance-covariance matrix of the random-effects}
\item{estimlink}{table containing the simulated values of each outcome and the corresponding estimated link function}
\item{epsY}{definite positive reals used to rescale the markers in (0,1) when the beta link function is used. By default, epsY=0.5.}
\item{linktype}{indicators of link function types: 0 for linear, 1 for beta, 2 for splines and 3 for thresholds}
\item{linknodes}{vector of nodes useful only for the 'splines' link functions}
%% idea0, idprob0,idg0,idcontr0,idcor0,Xnames2,na.action,pred_RE_Y,Ynames,nbnodes
}


\author{
Cecile Proust-Lima and Viviane Philipps

\email{cecile.proust-lima@inserm.fr}
}


\references{


Proust-Lima C, Philipps V, Liquet B (2015). Estimation of Extended Mixed Models Using Latent Classes and Latent Processes: the R package lcmm, Arxiv

Proust and Jacqmin-Gadda (2005). Estimation of linear mixed models with a mixture of distribution for the random-effects. Comput Methods Programs Biomed 78: 165-73.

Proust, Jacqmin-Gadda, Taylor, Ganiayre, and Commenges (2006). A
nonlinear model with latent process for cognitive evolution using multivariate longitudinal
data. Biometrics 62, 1014-24.

Proust-Lima, Dartigues and Jacqmin-Gadda (2011). Misuse of the linear mixed
model when evaluating risk factors of cognitive decline. Amer J Epidemiol 174(9): 1077-88.

Proust-Lima, Amieva, Jacqmin-Gadda (2013). Analysis of multivariate mixed longitudinal data: A flexible latent process approach. Br J Math Stat Psychol 66(3): 470-87.

Commenges, Proust-Lima, Samieri, Liquet (2012). A universal approximate cross-validation criterion and its
asymptotic distribution, Arxiv.
}


\seealso{

\code{\link{postprob}}, \code{\link{plot.multlcmm}}, \code{\link{predictL}}, \code{\link{predictY}} \code{\link{lcmm}}
}

\examples{
\dontrun{
data(data_Jointlcmm)
# Latent process mixed model for two curvilinear outcomes. Link functions are 
# aproximated by I-splines, the first one has 3 nodes (i.e. 1 internal node 8),
# the second one has 4 nodes (i.e. 2 internal nodes 12,25)

m1 <- multlcmm(Ydep1+Ydep2~1+Time*X2+contrast(X2),random=~1+Time,
subject="ID",randomY=TRUE,link=c("4-manual-splines","3-manual-splines"),
intnodes=c(8,12,25),data=data_Jointlcmm)

# to reduce the computation time, the same model is estimated using 
# a vector of initial values
m1 <- multlcmm(Ydep1+Ydep2~1+Time*X2+contrast(X2),random=~1+Time,
subject="ID",randomY=TRUE,link=c("4-manual-splines","3-manual-splines"),
intnodes=c(8,12,25),data=data_Jointlcmm, 
B=c(-1.071, -0.192,  0.106, -0.005, -0.193,  1.012,  0.870,  0.881,
  0.000,  0.000, -7.520,  1.401,  1.607 , 1.908,  1.431,  1.082,
 -7.528,  1.135 , 1.454 , 2.328, 1.052))


# output of the model
summary(m1)
# estimated link functions
plot(m1,which="linkfunction")
# variation percentages explained by linear mixed regression
VarExpl(m1,data.frame(Time=0))

#### Heterogeneous latent process mixed model with linear link functions 
#### and 2 latent classes of trajectory 
m2 <- multlcmm(Ydep1+Ydep2~1+Time*X2,random=~1+Time,subject="ID",
link="linear",ng=2,mixture=~1+Time,classmb=~1+X1,data=data_Jointlcmm,
B=c( 18,-20.77,1.16,-1.41,-1.39,-0.32,0.16,-0.26,1.69,1.12,1.1,10.8,
1.24,24.88,1.89))
# summary of the estimation
summary(m2)
# posterior classification
postprob(m2)
# longitudinal predictions in the outcomes scales for a given profile of covariates 
newdata <- data.frame(Time=seq(0,5,length=100),X1=rep(0,100),X2=rep(0,100),X3=rep(0,100))
predGH <- predictY(m2,newdata,var.time="Time",methInteg=0,nsim=20) 
head(predGH)
}
}

