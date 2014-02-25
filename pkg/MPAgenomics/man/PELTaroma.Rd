\name{PELTaroma}
\alias{PELTaroma}
\title{segmentation function}
\usage{
PELTaroma(dataSetName, normalTumorArray, chromosome = 1:22, Lambda = NULL,
  listOfFiles = NULL, onlySNP = TRUE, plot = TRUE, verbose = TRUE)
}
\arguments{
  \item{dataSetName}{name of the dataset (it must correpond
  to a folder name in rawData folder.).}

  \item{chromosome}{the chromosome which the CN signal must
  be extract.}

  \item{normalTumorArray}{only if you have normal and tumor
  profile in your data folder. A csv file or a data.frame
  with 2 columns: "normal" and "tumor". The first column
  contains the name of normal files and the second the
  names of associated tumor files.}

  \item{Lambda}{sequence of penalty parameters to test}

  \item{listOfFiles}{vector containing the names of the
  files in dataSetName to extract.}

  \item{onlySNP}{if TRUE only the copy-number at SNPs
  probes will be returned.}

  \item{plot}{if TRUE, plot the segmentation results}

  \item{verbose}{if TRUE print some informations}
}
\value{
a list containing \describe{ \item{copynumber}{A vector
containing the copynumber signal.} \item{segmented}{A
vector of the same size as copynumber containing the
segmented values.} \item{startPos}{The position of each
probes.} \item{chromosome}{A vector of the same size as
copynumber containing the chromosome number.}
\item{featureNames}{Names of the probes.}
\item{sampleNames}{The name of the signal.}
\item{segment}{A data.frame, each row is a different
segment with the chromosome, start position, end position,
number of probes in the signal and the value of the
segment.} }
}
\description{
This function launches the segmentation PELT (from package
changepoint) with a penalty value of lambda * log(n) for a
range of value for lambda. Then an optimal penalty value
will be choose in looking for a stabilization in the number
of segments according to the penalty values.
}
\author{
Quentin Grimonprez
}

