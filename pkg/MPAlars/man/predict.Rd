\name{predict}
\alias{predict}
\title{Prediction of response}
\usage{
  predict(x, Xnew, lambda, mode = "fraction")
}
\arguments{
  \item{x}{a LarsParth object}

  \item{Xnew}{a matrix (of size n*x@p) of covariates.}

  \item{lambda}{If mode =\"norm\", lambda represents the
  l1-norm of the coefficients with which we want to
  predict. If mode=\"fraction\", lambda represents the
  ratio (l1-norm of the coefficientswith which we want to
  predict)/(l1-norm maximal of the LarsPath object).}

  \item{mode}{"fraction" or "norm".}
}
\value{
  The predicted response
}
\description{
  Predict response of a new sample Xnew at a given level of
  penalty
}
\examples{
dataset=MPA.simul(50,10000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
result=MPA.lars(dataset$data[1:40,],dataset$response[1:40])
y=predict(result,dataset$data[41:50,],0.3,"fraction")
}
\author{
  Quentin Grimonprez
}

