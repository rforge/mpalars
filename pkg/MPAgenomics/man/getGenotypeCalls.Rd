\name{getGenotypeCalls}
\alias{getGenotypeCalls}
\title{Extract genotype calls from aroma files}
\usage{
getGenotypeCalls(dataSetName, chromosome, listOfFiles = NULL,
  verbose = TRUE)
}
\arguments{
  \item{dataSetName}{The name of the data-set folder (it
  must correpond to a folder name in rawData folder.)}

  \item{chromosome}{A vector containing the chromosomes for
  which the genotype call will be extracted.}

  \item{listOfFiles}{A vector containing the names of the
  files in dataSetName folder for which the genotype signal
  will be extracted (default is all the files).}

  \item{verbose}{If TRUE print some information
  (default=TRUE)}
}
\value{
a list of length the number of chromosome containing a
data.frame with columns: \describe{
\item{chromosome}{Chromosome of the signal.}
\item{position}{Positions associated with the genotype.}
\item{genotype}{Genotype calls corresponding to selected
files; the name of each column is the name of the
associated data file name.} \item{featureNames}{Names of
the probes.} }
}
\description{
Extract genotype calls from aroma files. It requires to
have executed the normalization process of aroma package.
}
\details{
You have to respect the aroma architecture. Your working
directory must contain rawData folder and callData folder.
To access easily the names of the files available in your
dataset, you can use the \link{getListOfFiles} function.
}
\examples{
#DO NOT EXECUTE
#fracB=getGenotypeCalls("data1",5)
}
\author{
Quentin Grimonprez
}

