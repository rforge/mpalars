\name{peltObject}
\alias{peltObject}
\title{create the list in parameters for \link{PELT} function}
\usage{
  peltObject(copynumber, chromosome, position,
    featureNames, sampleNames)
}
\arguments{
  \item{copynumber}{A vector containing the copy-number
  signal for one patient and one chromosome.}

  \item{chromosome}{Chromosome associated with the
  copy-number signal.}

  \item{position}{Position of the signal.}

  \item{featureNames}{Names of the probes (not necessary).}

  \item{sampleNames}{Name of the sample (not necessary).}
}
\value{
  a list in the right format for \link{PELT} function
}
\description{
  create the list in parameters for \link{PELT} function
}

