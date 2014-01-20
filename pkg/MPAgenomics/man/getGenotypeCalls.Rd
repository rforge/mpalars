\name{getGenotypeCalls}
\alias{getGenotypeCalls}
\title{Extract genotype calls from aroma files}
\usage{
  getGenotypeCalls(dataSetName, chromosome,
    listOfFiles = NULL, verbose = TRUE)
}
\arguments{
  \item{dataSetName}{Name of the dataset (it must correpond
  to a folder name in rawData folder.)}

  \item{chromosome}{A vector containing the chromosomes for
  which the genotype calls must be extract.}

  \item{listOfFiles}{A vector containing the names of the
  files in dataSetNames to extract.}

  \item{verbose}{If TRUE, print some informations.}
}
\value{
  a list of length the number of chromosome containing a
  data.frame with columns: \describe{
  \item{chromosome}{chromosome corresponding to the
  signal.} \item{position}{Positions associated to the
  allele B fraction.} \item{genotype}{Many columns (the
  number of selected files), the name of each column is the
  name of the associated data file name. It contains the
  genotype calls for a different profile.}
  \item{featureNames}{Names of the probes.} }
}
\description{
  Extract genotype calls from aroma files. It requires to
  have executed the normalization process of aroma package.
}
\details{
  You have to respect the aroma architecture. Your working
  directory must contain rawData folder and callData
  folder. To access easily the names of the files available
  in your dataset, you can use the \link{getListOfFiles}
  function.
}
\examples{
#DO NOT EXECUTE
#fracB=getGenotypeCalls("data1",5)
}

