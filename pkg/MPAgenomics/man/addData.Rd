\name{addData}
\alias{addData}
\title{Add a new data-set to the existing aroma architecture}
\usage{
addData(dataSetName, dataPath, chipType, verbose = TRUE)
}
\arguments{
  \item{dataSetName}{Name of the data-set folder to
  create.}

  \item{dataPath}{Path of the folder containing the data
  CEL files.}

  \item{chipType}{Name of the used chip.}

  \item{verbose}{Print additionnal information.}
}
\description{
Create a folder in "rawData" and copy the specified files
in this folder.
}
\author{
Quentin Grimonprez
}

