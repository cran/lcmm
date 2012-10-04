\name{lcmm-package}
\alias{lcmm-package}
\docType{package}
\title{
Estimation of Mixed Models, Latent Class Mixed Models and Joint Latent Class Mixed Models for different types of outcomes (quantitative, bounded, curvilinear and ordinal outcomes).
}
\description{
This package provides functions for the estimation of latent class mixed models (LCMM) and joint latent class mixed models (JLCM) using a maximum likelihood method. It also estimates mixed models for curvilinear and ordinal longitudinal outcomes (with or without latent classes of trajectories).
}
\details{
\tabular{ll}{
Package: \tab lcmm \cr
Type: \tab Package \cr
Version: \tab 1.5.7 \cr
Date: \tab 2012-10-04\cr
License: \tab GPL (>=2.0) \cr
LazyLoad: \tab yes \cr
}
The package includes for the moment the estimation of latent class mixed models for Gaussian longitudinal outcomes using \code{hlme} function, and for other quantitative, bounded quantitative (curvilinear) and discrete longitudinal outcomes using \code{lcmm} function, and joint latent class mixed models for a Gaussian longitudinal outcome and a right-censored (potentially left-truncated) time-to-event using \code{Jointlcmm} function. 

Please report to the maintainer any bug or comment regarding the package for future updates. 
}
\author{
Cecile Proust-Lima, Amadou Diakite and Benoit Liquet

Maintainer: <cecile.proust@isped.u-bordeaux2.fr>
}
\references{

Commenges, Liquet Proust-Lima (2012). Choice of prognostic estimators in joint models by estimating differences of
expected conditional {K}ullback-{L}eibler risks. Biometrics - in press

Lin, Turnbull, McCulloch and Slate (2002). Latent class models for
joint analysis of longitudinal biomarker and event process data: application to longitudinal
prostate-specific antigen readings and prostate cancer. Journal of the American Statistical
Association 97, 53-65.

Muthen and Shedden (1999). Finite mixture modeling with mixture outcomes using the EM algorithm. Biometrics 55, 463-9

Proust and Jacqmin-Gadda (2005). Estimation of linear mixed models with a mixture of distribution for the random-effects. Comput Methods Programs Biomed 78:165-73

Proust, Jacqmin-Gadda, Taylor, Ganiayre, and Commenges (2006). A
nonlinear model with latent process for cognitive evolution using multivariate longitudinal
data. Biometrics 62, 1014-24.

Proust-Lima, Dartigues and Jacqmin-Gadda (2011). Misuse of the linear mixed
model when evaluating risk factors of cognitive decline. Amer J Epidemiol 174(9), 1077-88

Proust-Lima and Taylor (2009). Development and validation of a dynamic prognostic
tool for prostate cancer recurrence using repeated measures of post-treatment PSA: a joint
modelling approach. Biostatistics 10, 535-49.

Proust-Lima, Sene, Taylor and Jacqmin-Gadda (2012). Joint latent class models of longitudinal and time-to-event data: a review. Statistical Methods in Medical Research - in press

Verbeke and Lesaffre (1996). A linear mixed-effects model with heterogeneity in
the random-effects population. Journal of the American Statistical Association 91, 217-21

}
\keyword{package}

