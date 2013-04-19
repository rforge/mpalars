\name{LarsPath-class}
\alias{LarsPath-class}
\title{Constructor of LarsPath class}
\description{
  This class stores the results of lars and fusion
  algorithms.
}
\details{
  \describe{ \item{nbStep}{Number of steps of the
  algorithm.} \item{variable}{List of vector of size
  "step+1". The i+1-th item contains the index of non-zero
  coefficients at the i-th step.} \item{coefficient}{List
  of vector of size "step+1". The i+1-th item contains the
  non-zero coefficients at the i-th step.}
  \item{lambda}{Vector of length "step+1", containing the
  L1-norm of the coefficients at each step.}
  \item{dropIndex}{Vector of length "step" containing the
  index of the dropped variable at the i-th step, 0 means
  no variable has been dropped at this step.}
  \item{addIndex}{Vector of length "step" containing the
  index of the added variable at the i-th step, 0 means no
  variable has been added at this step.}
  \item{mu}{Intercept.} \item{ignored}{A vector containing
  index of ignored variable during the algorithm.}
  \item{p}{Total number of covariates.} \item{fusion}{If
  TRUE, results from MPA.fusion function.} }
}

