\name{MPAlars-package}
\alias{MPAlars-package}
\alias{MPAlars}
\docType{package}
\title{
Lars algorithm
}

\description{
This package contains an implementation of the lars algorithm [1], for the lasso and fusion penalization. It works even if the number of covariate is greater than the number of individuals.

It has been designed for Multi Patients Analysis (MPA) of Genomic Markers.	

Its main purpose is the analysis of copy number signal but it can be used for other context.
}

\details{

\tabular{ll}{
Package: \tab MPAlars\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-05-10\cr
License: \tab GPL(>=2)\cr
}


}

\author{

Maintainer: Quentin Grimonprez <quentin.grimonprez@inria.fr>
}

\references{
[1] Efron, B., Hastie, T., Johnstone, I. and Tibshirani, R. 2004.  Least angle regression. Annals of Statistics, vol. 32 ,p407-499.
}

\keyword{ package }

\seealso{
}

\examples{
 dataset=MPA.simul(50,10000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2)) # simulate a data set
 result=MPA.lars(dataset$data,dataset$response) #launch lars algorithm
}
