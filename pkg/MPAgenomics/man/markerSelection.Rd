\name{markerSelection}
\alias{markerSelection}
\title{markers selection}
\usage{
markerSelection(dataSetName, dataResponse, chromosome = 1:22, signal = "CN",
  normalTumorArray, onlySNP = FALSE, nbFolds = min(length(dataResponse),
  10), loss = "logistic", plot = TRUE, ...)
}
\arguments{
  \item{dataSetName}{The name of the data-set folder.}

  \item{dataResponse}{A csv files or a data.frame with 2
  columns : "files" and "response". The column "files"
  contains the filename to extract and the second column
  the response associated to the file.}

  \item{chromosome}{A vector containing the number of the
  chromosomes for the SNPs selection.}

  \item{signal}{either \"CN\" or \"fracB\". corresponding
  to which signal will be analyzed (default=\”CN\”).}

  \item{normalTumorArray}{Only in the case of normal-tumor
  study. A csv file or a data.frame containing the mapping
  between normal and tumor files. The first column contains
  the name of normal files and the second the names of
  associated tumor files.}

  \item{onlySNP}{(only if signal=\"CN\"). If TRUE, only the
  SNPs probes are used (default=FALSE).}

  \item{nbFolds}{number of folds in the cross validation
  (default=10).}

  \item{loss}{either \"logistic\" (binary response) or
  \"linear\" (quantitative response), default is
  ”logistic”}

  \item{plot}{If TRUE, cross-validation mean squared error
  is plotted (default=TRUE).}

  \item{...}{Other parameters for HDlars function or glmnet
  function.}
}
\value{
a list containing length(chromosme) elements. Each element
is a list containing \describe{ \item{chr}{chromosome
corresponding to the signal.} \item{markers.index}{A vector
containing the index of all selected markers.}
\item{markers.position}{A vector containing the position of
all selected markers.} \item{markers.names}{A vector
containing the names of all selected markers.}
\item{coefficient}{A vector containing the coefficients of
all selected markers.} \item{intercept}{Intercept of the
model.} }
}
\description{
This function selects, for each chromosome, the most
relevant markers according to a response.
}
\details{
This function requires to use the aroma folder
architecture. In your working directory, there must have
the rawData folder and totalAndFracBData folder. This
function launches the lars algorithm on the CN or fracB
data and uses a cross-validation to select the most
appropriate solution.
}
\examples{
#DO NOT EXECUTE
# res=markerSelection("DataB",rnorm(27),chromosome=22,signal="CN",normalTumorArray,onlySNP=TRUE)
}
\author{
Quentin Grimonprez
}
\seealso{
HDPenReg
}

