\name{createArchitecture}
\alias{createArchitecture}
\title{Create aroma architecture and copy files}
\usage{
  createArchitecture(dataSetName, chipType, dataSetPath,
    chipFilesPath, path = ".", verbose = FALSE)
}
\arguments{
  \item{dataSetName}{name of the data set}

  \item{chipType}{type of the chip used for the data}

  \item{dataSetPath}{Path to the folder containing the data
  CEL files}

  \item{chipFilesPath}{Path to the folder containing the
  chip files}

  \item{path}{path where create folders}

  \item{verbose}{if TRUE, print details of the process}
}
\description{
  Create the architecture required by aroma.* packages and
  copy files into created folders.
}
\examples{
#createArchitecture("data1","GenomeWideSNP_6",dataSetPath="./data/",chipFilesPath="./chip/",path="./aroma",verbose=TRUE)
}
\seealso{
  copyChipFiles, copyDataFiles, createAromaArchitecture
}

