\name{addData}
\alias{addData}
\title{add a new data-set to the existing aroma architecture}
\usage{
addData(dataSetName, dataPath, chipType, verbose = TRUE)
}
\arguments{
  \item{dataSetName}{The name of the data-set folder to
  create.}

  \item{dataPath}{Path of the folder containing the data
  CEL files.}

  \item{chipType}{The name of the used chip.}
}
\description{
Create a folder in "rawData" and copy the specified files
in this folder.
}
\author{
Quentin Grimonprez
}

