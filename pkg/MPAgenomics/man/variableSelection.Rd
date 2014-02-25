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

  \item{plot}{if TRUE, plot the cross validation graphic.}

  \item{...}{Other parameters for HDlars function or glmnet
  function.}
}
\value{
a list containing \describe{ \item{variable}{A vector
containing the index of all selected variables.}
\item{coefficient}{A vector containing the coefficients of
all selected variables.} \item{intercept}{Intercept of the
model.} \item{fraction}{fraction of the l1 norm of
coefficient selected by cross validation.} }
}
\description{
This function selects the most relevant variables according
to a response.
}
\author{
Quentin Grimonprez
}

