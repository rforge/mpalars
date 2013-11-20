\name{HDPenReg-package}
\alias{HDPenReg-package}
\alias{HDPenReg}
\docType{package}
\title{
Algorithms for lasso and fused-lasso problems.
}

\description{
This package contains an implementation of the lars algorithm [1], for the lasso and fusion penalization. It works even if the number of covariate is greater than the number of individuals.
}

\details{

\tabular{ll}{
Package: \tab HDPenReg\cr
Type: \tab Package\cr
Version: \tab 0.88.1\cr
Date: \tab 2013-11-20\cr
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



\examples{
 dataset=simul(50,10000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2)) # simulate a data set
 result=HDlars(dataset$data,dataset$response) #launch lars algorithm
}
