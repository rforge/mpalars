\name{MPA.fusion}
\alias{MPA.fusion}
\title{Fusion algorithm}
\usage{
  MPA.fusion(X, y, maxSteps = 3 * min(dim(X)),
    eps = .Machine$double.eps^0.5, verbose = FALSE)
}
\arguments{
  \item{X}{the matrix (of size n*p) of the covariates.}

  \item{y}{a vector of length n with the response.}

  \item{maxSteps}{Maximal number of steps for lars
  algorithm.}

  \item{eps}{Tolerance of the algorithm.}

  \item{verbose}{If TRUE print details during the process.}
}
\value{
  An object of type LarsPath.
}
\description{
  fusion algorithm
}
\examples{
data=MPA.simul(50,10000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
result=MPA.fusion(data$data,data$response)
}
\author{
  Quentin Grimonprez
}

