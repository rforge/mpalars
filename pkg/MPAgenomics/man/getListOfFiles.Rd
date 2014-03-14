\name{getListOfFiles}
\alias{getListOfFiles}
\title{Get the contents of a data folder}
\usage{
getListOfFiles(dataSetName, chipType)
}
\arguments{
  \item{dataSetName}{The name of a data-set folder}

  \item{chipType}{The name of the used chip}
}
\value{
The filenames of all the files in
rawData/dataSetName/chipType
}
\description{
Get the cel files of the specified dataSetName
}
\details{


if chipType is not provided, the function returns the files
for the first chip (in the alphabetic order)
}
\author{
Quentin Grimonprez
}

