\name{Jointlcmm}
\alias{Jointlcmm}
\title{
Estimation of joint latent class models for longitudinal and time-to-event data
}
\description{
This function fits joint latent class mixed models for a longitudinal outcome and a right-censored (possibly left-truncated) time-to-event. In this first version, only the Gaussian longitudinal outcome is handled but other types of longitudinal outcomes should be handled soon (continuous -not Gaussian- as well as bounded quantitative, discrete and ordinal) similarly as in \code{lcmm} funtion. 
} 

\usage{
Jointlcmm(fixed, mixture, random, subject, classmb, ng=1, idiag=FALSE,
 nwg=FALSE, survival,hazard="Weibull",hazardtype="Specific", 
hazardnodes=NULL,TimeDepVar=NULL,data, B, convB=1e-4, convL=1e-4, 
convG=1e-4, maxiter=500, nsim=100, prior)
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
  \item{subject}{name of the covariate representing the grouping structure.
}
  \item{classmb}{
an optional one-sided formula describing the covariates in the class-membership multinomial logistic model. Covariates included are separated by \code{+}. No intercept should be included in this formula.
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
  \item{B}{
an optional vector containing the initial values for the parameters. The order in which the parameters are included is detailed in \code{details} section. 
If no vector is specified, a preliminary analysis involving the estimation of a standard linear mixed model is performed to choose initial values. 
Due to possible local maxima, the \code{B} vector should be specified and several different starting points should be tried.
}
\item{convB}{optional threshold for the convergence criterion based on the parameter stability. By default, convB=0.0001.
}
\item{convL}{optional threshold for the convergence criterion based on the log-likelihood stability. By default, convL=0.0001.
}
\item{convG}{optional threshold for the convergence criterion based on the derivatives. By default, convG=0.0001.
}

\item{maxiter}{optional maximum number of iterations for the Marquardt iterative algorithm. By default, maxiter=500.}

\item{nsim}{optional number of points for the predicted survival curves and predicted baseline risk curves. By default, nsim=100.
}
\item{prior}{name of the covariate containing the prior on the latent class membership. The covariate should be an integer with values in 0,1,...,ng. When there is no prior, the value should be 0. When there is a prior for the subject, the value should be the number of the latent class (in 1,...,ng).
}

\item{survival}{a two-sided formula object. The left side of the formula corresponds to a \code{surv()} object of type "counting" for right-censored and left-truncated data (example: \code{Surv(Time,EntryTime,Indicator)}) or of type "right" for right-censored data (example: \code{Surv(Time,Indicator)}). The right side of the formula includes names of covariates to include in the survival model with \code{mixture()} when the effect is class-specific (example: \code{Surv(Time,Indicator) ~ X1 + mixture(X2)} for a class-common effect of X1 and a class-specific effect of X2).}


\item{hazard}{optional family of hazard function assumed for the survival model. By default, "Weibull" specifies a Weibull baseline risk function. Other possibilities are "piecewise" for a piecewise constant risk function or "splines" for a cubic M-splines baseline risk function. For these two latter families, the number of nodes and the location of the nodes should be specified as well separated by \code{-}. The number of nodes is first entered followed by \code{-}, then the location is specified with "equi", "quant" or "manual" for respectively equidistant nodes, nodes at quantiles of the times of event distribution or interior nodes entered manually in argument \code{intnodes}. It is followed by \code{-} and finally "piecewise" or "splines" indicates the family of baseline risk function considered. Examples include "5-equi-splines" for M-splines with 5 equidistant nodes, "6-quant-piecewise" for piecewise constant risk over 5 intervals and nodes defined at the quantiles of the times of events distribution and "9-manual-splines" for M-splines risk function with 9 nodes, the vector of 7 interior nodes being entered in the argument \code{hazardnodes}.
 }


\item{hazardtype}{optional indicator for the type of baseline risk function when ng>1. By default "Specific" indicates a class-specific baseline risk function. Other possibilities are "PH" for a baseline risk function proportional in each latent class, and "Common" for a baseline risk function that is common over classes.}

\item{hazardnodes}{optional vector containing internal nodes if \code{splines} or \code{piecewise} is specified for the baseline hazard function in \code{hazard}.}


\item{TimeDepVar}{optional vector containing an intermediate time corresponding to a change in risk of event. This time-dependent covariate can only take the form of a time variable with the assumption that there is no effect before this time and a constant effect after on the risk of event (example: initiation of an additional treatment to account for).}


 \item{data}{
optional data frame containing the variables named in \code{fixed}, \code{mixture}, \code{random}, \code{classmb} and \code{subject}. 
}
}



\details{

A. BASELINE RISK FUNCTIONS

For the baseline risk functions, the following parameterizations were considered:

1. With the "Weibull" function: 2 parameters are necessary w_1 and w_2 so that the baseline risk function a_0(t) = w_1^2*w_2^2*(w_1^2*t)^(w_2^2-1). 

2. with the "piecewise" step function and nz nodes (y_1,...y_nz), nz-1 parameters are necesssary p_1,...p_nz-1 so that the baseline risk function a_0(t) = p_j^2 for y_j < t =< y_{j+1}.

3. with the "splines" function and nz nodes (y_1,...y_nz), nz+2 parameters are necessary s_1,...s_nz+2 so that the baseline risk function a_0(t) = sum_j s_j^2 M_j(t) where \{M_j\} is the basis of cubic M-splines.


B. THE VECTOR OF PARAMETERS B

The parameters in the vector of initial values \code{B} or in the vector of maximum likelihood estimates \code{best} are included in the following order: 
(1) ng-1 parameters are required for intercepts in the latent class membership model, and if covariates are included in \code{classmb}, ng-1 paramaters should be entered for each one; 
(2) parameters for the baseline risk function: 2 parameters for each Weibull, nz-1 for each piecewise constant risk and nz+2 for each splines risk; this number should be multiplied by ng if specific hazard is specified; otherwise, ng-1 additional proportional effects are expected if PH hazard is specified; otherwise nothing is added if common hazard is specified; 
(3) for all covariates in \code{survival}, ng parameters are required if the covariate is inside a \code{mixture()}, otherwise 1 parameter is required. Covariates parameters should be included in the same order as in \code{survival};
(4) for all covariates in \code{fixed}, one parameter is required if the covariate is not in \code{mixture}, ng paramaters are required if the covariate is also in \code{mixture}. Parameters should be included in the same order as in \code{fixed}; 
(5) the variance of each random-effect specified in \code{random} (including the intercept) 
if \code{idiag=TRUE} and the inferior triangular variance-covariance matrix of all the random-effects if \code{idiag=FALSE}; 
(6) only if \code{nwg=TRUE}, ng-1 parameters for class-specific proportional coefficients
 for the variance covariance matrix of the random-effects; 
(7) the variance of the residual error. 


C. CAUTION

Some caution should be made when using the program:


(1) As the log-likelihood of a latent class model can have multiple maxima, a careful choice of the initial values is crucial for ensuring convergence toward the global maximum. 
The program can be run without entering the vector of initial values (see point 2). 
However, we recommend to systematically enter initial values in \code{B} and try different sets of initial values. 


(2) The automatic choice of initial values that we provide requires the estimation of a preliminary linear mixed model. The user should be aware that first, this preliminary analysis can take time for large datatsets and second, 
that the generated initial values can be very not likely and even may converge slowly to a local maximum. 
This is a reason why specification of initial values in \code{B} should be preferred.


(4) Convergence criteria are very strict as they are based on derivatives of the log-likelihood in addition to the parameter and log-likelihood stability. 
In some cases, the program may not converge and reach the maximum number of iterations fixed at 70. 
In this case, the user should check that parameter estimates at the last iteration are not on the boundaries of the parameter space. 
If the parameters are on the boundaries of the parameter space, the identifiability of the model should be assessed. 
If not, the program should be run again with other initial values.
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
\item{pprob}{table of posterior classification and posterior individual class-membership probabilities based on the longitudinal data and the time-to-event data}
\item{pprobY}{table of posterior classification and posterior individual class-membership probabilities based only on the longitudinal data}
\item{Xnames}{list of covariates included in the model - for use in function \code{\link{plot.predict.hlme}}}
\item{predRE}{table containing individual predictions of the random-effects: a column per random-effect, a line per subject}
\item{cholesky}{vector containing the estimates of the Cholesky transformed parameters of the variance-covariance matrix of the random-effects}
\item{CIstat}{Statistic of the Score Test for the conditional independence assumption of the longitudinal and survival data given the latent class structure. Under the null hypothesis, the statistics is a Chi-square with p degrees of freedom where p indicates the number of random-effects in the longitudinal mixed model. See Jacqmin-Gadda and Proust-Lima (2009) for more details.}
\item{predSurv}{table of predictions giving for the window of times to event (called "time"), the predicted baseline risk function in each latent class (called "RiskFct") and the predicted cumulative baseline risk function in each latent class (called "CumRiskFct").}
}

\author{
Cecile Proust Lima and Amadou Diakite 

Maintainer: <cecile.proust@isped.u-bordeaux2.fr>
}

\references{
Lin, H., Turnbull, B. W., McCulloch, C. E. and Slate, E. H. (2002). Latent class models for
joint analysis of longitudinal biomarker and event process data: application to longitudinal
prostate-specific antigen readings and prostate cancer. Journal of the American Statistical
Association 97, 53-65.

Proust-Lima, C. and Taylor, J. (2009). Development and validation of a dynamic prognostic
tool for prostate cancer recurrence using repeated measures of post-treatment PSA: a joint
modelling approach. Biostatistics 10, 535-49.


Jacqmin-Gadda, H. and Proust-Lima, C. (2010). Score test for conditional independence between longitudinal outcome and time-to-event given the classes in the joint latent class model. Biometrics 66(1), 11-9



}
\seealso{
\code{\link{postprob}}, \code{\link{plot.survival}}, \code{\link{plot.baselinerisk}}, \code{\link{plot.predict}}
}
\examples{

\dontrun{

##### Example of a joint latent class model estimated for a varying number
# of latent classes: 
# The linear mixed model includes a subject (ID) and class-specific linear 
# trend (intercept and Time in fixed, random and mixture components) and
# a common effect of X1 and its interaction with time (X1_Time) over classes 
# (in fixed).
# The variance of the random intercept and slopes are assumed to be equal 
# over classes (nwg=F).
# The covariate X3 predicts the  class membership (in classmb). 
# The baseline hazard function is modelled with cubic M-splines -5 
# equidistant nodes- (in hazard) and a proportional hazard over classes is 
# assumed (in hazardtype). Covariates X1 and X2 predict the risk of event 
# (in survival) with a common effect over classes for X1 and a class-specific
# effect of X2.
# !CAUTION: for illustration, only default initial values where used but 
# other sets of initial values should be tried to ensure convergence towards
# the global maximum.



#### data loading
data(data_Jointlcmm)

#### estimation with 1 latent class (ng=1): independent models for the 
# longitudinal outcome and the time of event
m1 <- Jointlcmm(fixed= Ydep1~Time+X1+X1_Time,random=~Time,subject='ID'
,survival = Surv(Tevent,Event)~ X1+X2 ,hazard="5-equi-splines"
,hazardtype="PH",ng=1,data=data_Jointlcmm)
#Goodness-of-fit statistics for m1:
#     maximum log-likelihood: -3938.91; AIC: 7911.82; BIC: 7974.79

#### estimation with 2 latent classes (ng=2)
m2 <- Jointlcmm(fixed= Ydep1~Time+X1+X1_Time,mixture=~Time,random=~Time,
classmb=~X3,subject='ID',nwg=F,survival = Surv(Tevent,Event)~X1+mixture(X2)
,hazard="5-equi-splines",hazardtype="PH",ng=2,data=data_Jointlcmm)
summary(m2)
#Goodness-of-fit statistics for m2:
#     maximum log-likelihood: -3913.96;  AIC: 7873.93; BIC: 7959.12

#### estimation with 3 latent classes (ng=3)
m3 <- Jointlcmm(fixed= Ydep1~Time+X1+X1_Time,mixture=~Time,random=~Time,
classmb=~X3,subject='ID',survival = Surv(Tevent,Event)~ mixture(X1)+X2
,hazard="5-equi-splines",hazardtype="PH",ng=3,data=data_Jointlcmm)
summary(m3)
#Goodness-of-fit statistics for m3:
#     maximum log-likelihood: -3885.94 ;  AIC: 7829.89; BIC: 7937.29

#### estimation with 4 latent classes (ng=4)
m4 <- Jointlcmm(fixed= Ydep1~Time+X1+X1_Time,mixture=~Time,random=~Time,
classmb=~X3,subject='ID',survival = Surv(Tevent,Event)~mixture(X1)+X2,
hazard="5-equi-splines",hazardtype="PH",ng=4,data=data_Jointlcmm)
summary(m4)
#Goodness-of-fit statistics for m4:
#     maximum log-likelihood: -3882.38; AIC: 7834.75; BIC: 7964.38

##### The model with 3 latent classes is retained according to the BIC  
##### and the conditional independence assumption is not rejected at
##### the 5% level. 
# posterior classification
postprob(m3)
# Class-specific predicted baseline risk & survival functions in the 
# 3-class model retained (for the reference value of the covariates) 
plot.baselinerisk(m3)
plot.baselinerisk(m3,ylim=c(0,5))
plot.survival(m3)
# class-specific predicted trajectories in the 3-class model retained 
# (with characteristics of subject ID=193)
data <- data_Jointlcmm[data_Jointlcmm$ID==193,]
plot.predict(m3,var.time="Time",newdata=data)
############## end of example ##############

}
}

  
