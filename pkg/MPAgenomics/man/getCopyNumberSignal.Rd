\name{getCopyNumberSignal}
\alias{getCopyNumberSignal}
\title{Extract copy-number signal from aroma files}
\usage{
  getCopyNumberSignal(dataSetName, chromosome,
    normalTumorArray, onlySNP = FALSE, listOfFiles = NULL,
    verbose = TRUE)
}
\arguments{
  \item{dataSetName}{name of the dataset (it must correpond
  to a folder name in rawData folder.).}

  \item{chromosome}{A vector containing the chromosomes for
  which the CN signal must be extract.}

  \item{normalTumorArray}{only if you have normal and tumor
  profile in your data folder. A csv file or a data.frame
  with 2 columns: "normal" and "tumor". The first column
  contains the name of normal files and the second the
  names of associated tumor files.}

  \item{onlySNP}{if TRUE only the copy-number at SNPs
  probes will be returned.}

  \item{listOfFiles}{vector containing the names of the
  files in dataSetName to extract.}

  \item{verbose}{if TRUE print some informations.}
}
\value{
  a list of length the number of chromosome containing a
  data.frame with columns: \describe{
  \item{chromosome}{chromosome corresponding to the
  signal.} \item{position}{Positions associated to the
  allele B fraction.} \item{copynumber}{Many columns (the
  number of selected files), the name of each column is the
  name of the associated data file name. It contains the CN
  signal for a different profile.}
  \item{featureNames}{Names of the probes.} }
}
\description{
  Extract copy-number signals from aroma files. It requires
  to have executed the normalization process of aroma
  package.
}
\details{
  You have to respect the aroma architecture. Your working
  directory must contain rawData folder and
  totalAndFracBData folder. To access easily the names of
  the files available in your dataset, you can use the
  \link{getListOfFiles} function.
}
\examples{
#DO NOT EXECUTE
#C=getCopyNumberSignal("data1",5,normalTumorArray,TRUE)
#C=getCopyNumberSignal("data2",5,onlySNP=TRUE)
}

