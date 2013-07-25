\name{EMlasso}
\alias{EMlasso}
\title{EM algorithm for lasso penalty}
\usage{
  EMlasso(X, y, lambda, maxSteps = 3 * max(dim(X)),
    eps = .Machine$double.eps^0.5)
}
\arguments{
  \item{X}{the matrix (of size n*p) of the covariates.}

  \item{y}{a vector of length n with the response.}

  \item{lambda}{a sequence of l1 penalty regularization
  term. If no sequence is provided, the function computes
  his own sequence.}

  \item{maxSteps}{Maximal number of steps for lars
  algorithm.}

  \item{eps}{Tolerance of the algorithm.}
}
\value{
  An object of type LarsPath.
}
\description{
  EM algorithm for lasso penalty
}
\examples{
dataset=simul(50,10000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
result=EMlasso(dataset$data,dataset$response)
}
\author{
  Quentin Grimonprez
}

