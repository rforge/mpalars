\name{getFracBSignal}
\alias{getFracBSignal}
\title{Extract allele B fraction signal from aroma files}
\usage{
getFracBSignal(dataSetName, chromosome, normalTumorArray, listOfFiles = NULL,
  verbose = TRUE)
}
\arguments{
  \item{dataSetName}{Name of the dataset (it must correpond
  to a folder name in rawData folder.)}

  \item{chromosome}{A vector containing the chromosomes for
  which the allele B fraction signal must be extract.}

  \item{normalTumorArray}{Only if you have normal and tumor
  profile in your data folder. A csv file or a data.frame
  with 2 columns: "normal" and "tumor". The first column
  contains the name of normal files and the second the
  names of associated tumor files.}

  \item{listOfFiles}{A vector containing the names of the
  files in dataSetNames to extract.}

  \item{verbose}{If TRUE, print some informations.}
}
\value{
a list of length the number of chromosome containing a list
of two elements (normal and tumor) containing a data.frame
with columns: \describe{ \item{chromosome}{chromosome
corresponding to the signal.} \item{position}{Positions
associated to the allele B fraction.} \item{fracB}{Many
columns (the number of selected files), the name of each
column is the name of the associated data file name. It
contains the fracB signal for a different profile.}
\item{featureNames}{Names of the probes.} }
}
\description{
Extract allele B fraction signals from aroma files. It
requires to have executed the normalization process of
aroma package.
}
\details{
You have to respect the aroma architecture. Your working
directory must contain rawData folder and totalAndFracBData
folder. In case of single array study only the fracB$tumor
matrix is returned. To access easily the names of the files
available in your dataset, you can use the
\link{getListOfFiles} function.
}
\examples{
#DO NOT EXECUTE
#fracB=getFracBSignal("data1",5,normalTumorArray)
#fracB=getFracBSignal("data2",5)
}
\author{
Quentin Grimonprez
}

