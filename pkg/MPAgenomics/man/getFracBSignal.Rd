\name{getFracBSignal}
\alias{getFracBSignal}
\title{Extract allele B fraction signal from aroma files}
\usage{
getFracBSignal(dataSetName, chromosome, normalTumorArray, listOfFiles = NULL,
  verbose = TRUE)
}
\arguments{
  \item{dataSetName}{The name of the data-set folder (it
  must correspond to a folder name in rawData folder.)}

  \item{chromosome}{A vector containing the chromosomes for
  which the allele B fraction signal must be extract.}

  \item{normalTumorArray}{Only in the case of normal-tumor
  study. A csv file or a data.frame containing the mapping
  between normal and tumor files The first column contains
  the name of normal files and the second the names of
  associated tumor files.}

  \item{listOfFiles}{A vector containing the names of the
  files in dataSetName folder for which the allele B
  fraction profiles will be extracted (default is all the
  files).}

  \item{verbose}{If TRUE print some information
  (default=TRUE).}
}
\value{
a list of length the number of chromosomes containing a
list of two elements (normal and tumor) containing a
data.frame with columns: \describe{
\item{chromosome}{Chromosome of the signal.}
\item{position}{Positions associated with the allele B
fraction.} \item{fracB}{Allele B fraction profiles of
selected files; the name of each column is the name of the
associated data file name.} \item{featureNames}{Names of
the probes.} }
}
\description{
Extract allele B fraction signals from aroma files. It
requires to have executed the normalization process
suggested by aroma packages, by using
\link{signalPreProcess} for example.
}
\details{
The aroma architecture must be respected. The working
directory must contain rawData folder and totalAndFracBData
folder. To easily access the names of the files available
in a dataset, one can use the \link{getListOfFiles}
function.
}
\examples{
#DO NOT EXECUTE before reading the vignette
#fracB=getFracBSignal("data1",5,normalTumorArray)
#fracB=getFracBSignal("data2",5)
}
\author{
Quentin Grimonprez
}

