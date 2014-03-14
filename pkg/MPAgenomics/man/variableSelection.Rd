\name{variableSelection}
\alias{variableSelection}
\title{SNPs selection}
\usage{
variableSelection(dataMatrix, dataResponse,
  nbFolds = min(length(dataResponse), 10), loss = "logistic", plot = TRUE,
  ...)
}
\arguments{
  \item{dataMatrix}{Matrix containing the data, each row is
  a different sample.}

  \item{dataResponse}{response associated to the data.}

  \item{nbFolds}{number of folds in the cross validation.}

  \item{loss}{either \"logistic\" (binary response) or
  \"linear\" (quantitative response).}

  \item{plot}{If TRUE plot cross-validation mean squared
  error (default=TRUE).}

  \item{...}{spplementary arguments for cv.glmnet function
  in case of logistic loss or for HDlars function for
  linear loss.}
}
\value{
a list containing \describe{ \item{variable}{A vector
containing the index of all selected variables.}
\item{coefficient}{A vector containing the coefficients of
all selected variables.} \item{intercept}{Intercept of the
model.} }
}
\description{
This function selects the most relevant variables according
to a response.
}
\author{
Quentin Grimonprez
}

