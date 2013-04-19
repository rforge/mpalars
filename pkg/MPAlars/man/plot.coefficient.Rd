\name{plot.coefficient}
\alias{plot.coefficient}
\title{Plot of coefficients}
\usage{
  plot.coefficient(x, step, ylab = "coefficients",
    xlab = "variables", ...)
}
\arguments{
  \item{x}{A LarsPath object.}

  \item{step}{The step at which you want to plot the
  coefficients.}

  \item{ylab}{Name of the y axis.}

  \item{xlab}{Name of the x axis.}

  \item{...}{Other plot arguments.}
}
\description{
  Plot of the coefficients of a step
}
\examples{
data=MPA.simul(50,10000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
result=MPA.fusion(data$data,data$response)
plot.coefficient(result,45) #plot coefficients at the step 45
}

