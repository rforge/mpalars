\name{createArchitecture}
\alias{createArchitecture}
\title{Create aroma architecture and copy files}
\usage{
createArchitecture(dataSetName, chipType, dataSetPath, chipFilesPath,
  path = ".", verbose = FALSE, tags = NULL)
}
\arguments{
  \item{dataSetName}{The name of the data-set folder to
  create}

  \item{chipType}{The name of the used chip}

  \item{dataSetPath}{Path to the folder containing the data
  CEL files}

  \item{chipFilesPath}{Path to the folder containing the
  chip files}

  \item{path}{Path where the architecture should be created
  (default=”.”)}

  \item{verbose}{Print information during the process
  (default=FALSE)}

  \item{Common}{tag which appears in the different file
  names (cdf, ugp, ufl) of the chip. For no tag, use
  tags=NULL (default = NULL). See details for more
  information.}
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
createArchitecture("datatest1","GenomeWideSNP_6",dataSetPath="./celPATH",chipFilesPath="./chipPATH",path=".",verbose=TRUE,tags="Full")
}
\author{
Quentin Grimonprez
}
\seealso{
copyChipFiles, copyDataFiles, createAromaArchitecture
}

