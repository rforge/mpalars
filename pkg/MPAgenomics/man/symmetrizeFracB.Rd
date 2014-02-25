\name{symmetrizeFracB}
\alias{symmetrizeFracB}
\title{symmetrize an allele B fraction signal}
\usage{
symmetrizeFracB(fracB)
}
\arguments{
  \item{fracB}{a vector containing an allele B fraction
  signal.}
}
\value{
a vector containing the symmetrized signal.
}
\description{
The fraction allele B signal is the ratio between the
signal from the B allele and the total signal. The
symmetrization of the fraction allele B signal x is :
2*abs(x-0.5).
}
\examples{
signalA=abs(rnorm(100))
signalB=abs(rnorm(100))
signalFracB=signalA/(signalA+signalB)

symFracB=symmetrizeFracB(signalFracB)
}
\author{
Quentin Grimonprez
}

