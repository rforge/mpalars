\name{EMfusedlasso}
\alias{EMfusedlasso}
\title{EM algorithm for fused lasso penalty}
\usage{
EMfusedlasso(X, y, lambda1, lambda2, maxSteps = 10000, burn = 50,
  intercept = TRUE, model = "linear", eps = 1e-05, eps0 = 1e-12,
  epsCG = 1e-08)
}
\arguments{
  \item{X}{the matrix (of size n*p) of the covariates.}

  \item{y}{a vector of length n with the response.}

  \item{lambda1}{a positive real. Parameter associated with
  the lasso penalty.}

  \item{lambda2}{a positive real. Parameter associated with
  the fusion penalty.}

  \item{maxSteps}{Maximal number of steps for EM
  algorithm.}

  \item{burn}{Number of steps before regrouping some
  variables in segment.}

  \item{model}{"linear" or "logistic"}

  \item{intercept}{If TRUE, there is an intercept in the
  model.}

  \item{eps}{tolerance for convergence of the EM
  algorithm.}

  \item{eps0}{epsilon to add to avoid 0 to denominator of
  fraction.}

  \item{epsCG}{tolerance for convergence of the conjugate
  gradient.}
}
\value{
A list containing : \describe{ \item{nbStep}{Vector
containing the number of steps of the algorithm for every
lambda.} \item{variable}{List of vector of size "step+1".
The i+1-th item contains the index of non-zero coefficients
at the i-th step.} \item{coefficient}{List of vector of
size "step+1". The i+1-th item contains the non-zero
coefficients at the i-th step.} \item{lambda}{Vector of
length "step+1", containing the lambda at each step.}
\item{mu}{Intercept.} }
}
\description{
EM algorithm for fused lasso penalty
}
\examples{
dataset=simul(50,100,0.4,1,10,matrix(c(0.1,0.9,0.02,0.02),nrow=2))
result=EMfusedlasso(dataset$data,dataset$response,1,1)
}
\author{
Quentin Grimonprez
}

