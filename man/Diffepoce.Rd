\name{Diffepoce}
\alias{Diffepoce}
\title{
Computation of the difference of EPOCE estimators and its 95% tracking interval between two joint latent class models estimated using \code{Jointlcmm}
}
\description{
This functions computes the difference of CVPOL or MPOL and its 95% tracking interval between two joint latent class models estimated using \code{Jointlcmm} and evaluated using \code{epoce}. Difference in CVPOL is computed when EPOCE was previously estimated on the same dataset as used for estimation (using an approximated cross-validation), and difference in MPOL is computed when EPOCE was previously estimated on an external dataset. 
}
\usage{
Diffepoce(epoceM1,epoceM2)
}
\arguments{
  \item{epoceM1}{a first object inheriting from class \code{epoce}
  }
  \item{epoceM2}{a second object inheriting from class \code{epoce}
  }
}
\details{
From the EPOCE estimates and the individual contributions to the prognostic observed log-likelihood obtained on the same dataset with \code{epoce} for two different joint latent class models, the difference of CVPOL (or MPOL) and its 95% tracking interval is computed. The 95% tracking interval is: 

Delta(MPOL) +/- qnorm(0.95)*sqrt(VARIANCE) for an external dataset
Delta(CVPOL) +/- qnorm(0.95)*sqrt(VARIANCE) for the dataset used in \code{Jointlcmm}

where Delta(CVPOL) (Delta(MPOL)) is the difference of CVPOL or MPOL of the two joint latent class models, and VARIANCE is the empirical variance of the difference of individual contributions to the prognostic observed log-likelihoods of the two joint latent class models.

See Commenges et al. (2012) for further details.
}
\value{
\item{call.Jointlcmm1}{the \code{Jointlcmm} call for the first \code{epoce} object in epoceM1
}
\item{call.Jointlcmm2}{the \code{Jointlcmm} call for the second \code{epoce} object in epoceM2
}
\item{call}{the matched call}
\item{DiffEPOCE}{Dataframe containing, for each prediction time s, the difference in either MPOL or CVPOL depending on the dataset used, and the 95% tracking bands (TIinf and TIsup)
}
\item{new.data}{a boolean for internal use only, which is FALSE if computation is done on the same data as for \code{Jointlcmm} estimation, and TRUE otherwise.
}
}
\references{
Commenges, Liquet and Proust-Lima (2012). Choice of prognostic estimators in joint models by estimating differences of
expected conditional {K}ullback-{L}eibler risks. Biometrics - in press

Proust-Lima, Sene, Taylor and Jacqmin-Gadda (2012). Joint latent class models of longitudinal and time-to-event data: a review. Statistical Methods in Medical Research - in press
}
\author{
Cecile Proust-Lima and Amadou Diakite
}

\seealso{
\code{\link{Jointlcmm}},\code{\link{epoce}},\code{\link{summary.Diffepoce}}
}
\examples{
#### estimation with 2 latent classes (ng=2)
data(data_Jointlcmm)
m2 <- Jointlcmm(fixed= Ydep1~Time*X1,random=~Time,mixture=~Time,subject='ID'
,survival = Surv(Tevent,Event)~ X1+X2 ,hazard="Weibull"
,hazardtype="PH",ng=2,data=data_Jointlcmm,
B=c( 0.7608, -9.4974,  1.0242,  1.4331,  0.1063 , 0.6714, 10.4679, 11.3178,
 -2.5671, -0.5386,  1.4616, -0.0605,  0.9489,  0.1020,  0.2079,  1.5045))
m1 <- Jointlcmm(fixed= Ydep1~Time*X1,random=~Time,subject='ID'
,survival = Surv(Tevent,Event)~ X1+X2 ,hazard="Weibull"
,hazardtype="PH",ng=1,data=data_Jointlcmm,
B=c(-7.6634,  0.9136,  0.1002,  0.6641, 10.5675, -1.6589,  1.4767, -0.0806,
  0.9240,0.5643,  1.2277,  1.5004))

## EPOCE computation for predictions times from 1 to 6 on the dataset used for estimation of m.
VecTime <- c(1,3,5,7,9,11,13,15)
cvpol1 <- epoce(m1,var.time="Time",pred.times=VecTime)
cvpol1
cvpol2 <- epoce(m2,var.time="Time",pred.times=VecTime)
cvpol2
DeltaEPOCE <- Diffepoce(cvpol1,cvpol2)
summary(DeltaEPOCE)
plot(DeltaEPOCE,bty="l")
}

