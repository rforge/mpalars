\name{getListOfFiles}
\alias{getListOfFiles}
\title{Get the contents of a data folder}
\usage{
  getListOfFiles(dataSetName, chipType)
}
\arguments{
  \item{dataSetName}{Name of a folder in rawData folder}

  \item{chipType}{Name of the chip type used for the data}
}
\value{
  The filenames of all the files in
  rawData/dataSetName/chipType
}
\description{
  Get the contents of a data folder
}
\details{
  if chipType is not provided, the function returns the
  files for the first chip (in the alphabetic order)
}

