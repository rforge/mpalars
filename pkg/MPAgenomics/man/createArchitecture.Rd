\name{createArchitecture}
\alias{createArchitecture}
\title{Create aroma architecture and copy files}
\usage{
createArchitecture(dataSetName, chipType, dataSetPath, chipFilesPath,
  path = ".", verbose = FALSE, tags = NULL)
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

  \item{tags}{appear in the different file name (cdf, ugp,
  ufl) of the chip. For no tag, use tags=NULL (default =
  NULL). See details for more information.}
}
\description{
Create the architecture required by aroma.* packages and
copy files into created folders.
}
\details{
All the cdf chip file names must follow the following rule
: <chipType>,<Tags>.cdf Multiples tags must be separated by
a comma. If there is no tags, the pattern is <chipType>.cdf
}
\examples{
#createArchitecture("data1","GenomeWideSNP_6",dataSetPath="./data/",chipFilesPath="./chip/",path="./aroma",verbose=TRUE)
}
\author{
Quentin Grimonprez
}
\seealso{
copyChipFiles, copyDataFiles, createAromaArchitecture
}

