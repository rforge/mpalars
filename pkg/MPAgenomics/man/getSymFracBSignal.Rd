\name{getSymFracBSignal}
\alias{getSymFracBSignal}
\title{Extract symmetrized allele B fraction signal from aroma files}
\usage{
getSymFracBSignal(dataSetName, file, chromosome, normalTumorArray,
  verbose = TRUE)
}
\arguments{
  \item{dataSetName}{The name of the data-set folder (it
  must correpond to a folder name in rawData folder.)}

  \item{chromosome}{A vector with the chromosomes for which
  the symetrized signal will be extracted.}

  \item{normalTumorArray}{Only if you have normal and tumor
  profile in your data folder. A csv file or a data.frame
  with 2 columns: "normal" and "tumor". The first column
  contains the name of normal files and the second the
  names of associated tumor files.}

  \item{file}{The name of the file in dataSetName to
  extract.}

  \item{verbose}{If TRUE, print some informations.}
}
\value{
a list of length the number of chromosome containing a
data.frame with columns: \describe{
\item{chromosome}{chromosome corresponding to the signal.}
\item{position}{Positions associated to the allele B
fraction.} \item{fracB}{One column named by the data file
name. It contains the symmetrized allele B fraction signal
for the specified profile.} \item{featureNames}{Names of
the probes.} }
}
\description{
Extract symmetrized allele B fraction signals from aroma
files. It requires to have executed the normalization
process of aroma package.
}
\details{
You have to respect the aroma architecture. Your working
directory must contain rawData folder and totalAndFracBData
folder. In case of single array study only the fracB$tumor
matrix is returned.
}
\examples{
#DO NOT EXECUTE
#fracB=getSymFracBSignal("data1",5,normalTumorArray)
#fracB=getSymFracBSignal("data2",5)
}
\author{
Quentin Grimonprez
}

