\name{markerSelection}
\alias{markerSelection}
\title{markers selection}
\usage{
  markerSelection(dataSetName, dataResponse,
    chromosome = 1:22, signal = "CN", normalTumorArray,
    onlySNP = FALSE,
    nbFolds = min(length(dataResponse), 10),
    loss = "logistic", plot = TRUE, ...)
}
\arguments{
  \item{dataSetName}{Name of the dataset folder in rawData
  .}

  \item{dataResponse}{A csv files or a data.frame with 2
  columns : "files" and "response". The column "files"
  contains the filename to extract and the second column
  the response associated to the file.}

  \item{chromosome}{chromosome used in the study.}

  \item{signal}{either \"CN\" or \"fracB\". The type of
  signal to import.}

  \item{normalTumorArray}{only if you have normal and tumor
  profile in your data folder. A csv file or a data.frame
  with 2 columns: "normal" and "tumor". The first column
  contains the name of normal files and the second the
  names of associated tumor files.}

  \item{onlySNP}{(only if signal=\"CN\"). If TRUE, only the
  SNP markers will be used.}

  \item{nbFolds}{number of folds in the cross validation.}

  \item{loss}{either \"logistic\" (binary response) or
  \"linear\" (quantitative response).}

  \item{plot}{if TRUE, plot the cross validation graphic.}

  \item{...}{Other parameters for HDlars function or glmnet
  function.}
}
\value{
  a list containing length(chromosme) elements. Each
  element is a list containing \describe{
  \item{chr}{chromosome corresponding to the signal.}
  \item{markers.index}{A vector containing the index of all
  selected markers.} \item{markers.position}{A vector
  containing the position of all selected markers.}
  \item{markers.names}{A vector containing the names of all
  selected markers.} \item{coefficient}{A vector containing
  the coefficients of all selected markers.}
  \item{intercept}{Intercept of the model.} }
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
\seealso{
  HDPenReg
}

