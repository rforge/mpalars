\name{coeff}
\alias{coeff}
\title{get coefficients at a given step.}
\usage{
  coeff(x, step)
}
\arguments{
  \item{x}{A LarsPath object.}

  \item{step}{The step at which you want to get the
  coefficients.}
}
\value{
  a vector of size p containing the value of coefficients
  at the desired step.
}
\description{
  Get the vector of coefficients at a given step
}
\examples{
dataset=simul(50,10000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
result=HDfusion(dataset$data,dataset$response)
coeff(result,result@nbStep) #get the coefficients
}

