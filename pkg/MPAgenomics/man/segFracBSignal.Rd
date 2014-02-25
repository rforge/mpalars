\name{segFracBSignal}
\alias{segFracBSignal}
\title{segmentation function for the allele B fraction}
\usage{
segFracBSignal(dataSetName, normalTumorArray, chromosome = 1:22,
  Lambda = NULL, listOfFiles = NULL, savePlot = TRUE, verbose = TRUE)
}
\arguments{
  \item{dataSetName}{name of the dataset (it must correpond
  to a folder name in rawData folder.).}

  \item{normalTumorArray}{only if you have normal and tumor
  profile in your data folder. A csv file or a data.frame
  with 2 columns: "normal" and "tumor". The first column
  contains the name of normal files and the second the
  names of associated tumor files.}

  \item{chromosome}{A vector containing the chromosome
  which the CN signal must be extract.}

  \item{Lambda}{sequence of penalty parameters to test}

  \item{listOfFiles}{vector containing the names of the
  files in dataSetName to extract.}

  \item{savePlot}{if TRUE, plot the segmentation results}

  \item{verbose}{if TRUE print some informations}
}
\value{
a data.frame where each row correspond to a different
segment with columns : \describe{ \item{sampleNames}{The
name of the signal.} \item{chromosome}{A vector of the same
size as copynumber containing the chromosome number.}
\item{chromStart}{The starting position of a segment.}
\item{chromEnd}{The ending position of a segment.}
\item{probes}{The number of probes in the segment.}
\item{means}{Means of the segment.} }
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

