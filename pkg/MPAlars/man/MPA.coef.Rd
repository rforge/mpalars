\name{MPA.coef}
\alias{MPA.coef}
\title{get coefficients at a given step.}
\usage{
  MPA.coef(x, step)
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
dataset=MPA.simul(50,10000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
result=MPA.fusion(dataset$data,dataset$response)
MPA.coef(result,45) #get the coefficients
}

