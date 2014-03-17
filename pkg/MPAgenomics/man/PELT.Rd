\name{PELT}
\alias{PELT}
\title{segmentation function}
\usage{
PELT(signal, Lambda, position = NULL, plot = TRUE, verbose = TRUE)
}
\arguments{
  \item{signal}{a vector containing the signal.}

  \item{Lambda}{A vector containing all the penalization
  values to test for the segmentation. If no values are
  provided, default values will be used}

  \item{position}{A vector containing the position of all
  elements of the signal (not necessary)}

  \item{plot}{if TRUE, plot the segmentation results}

  \item{verbose}{if TRUE print some informations}
}
\value{
a list containing \describe{ \item{signal}{A vector
containing the signal.} \item{segmented}{A vector of the
same size as signal containing the segmented values.}
\item{startPos}{The position of each probe.}
\item{segment}{A data.frame that summarizes the results of
the segmentation. Each row is a different segment with the
start position, end position, number of points in the
signal and the value of the segment.} }
}
\description{
This function launches the segmentation PELT (from package
changepoint) with a penalty value of lambda * log(n) for a
range of value for lambda. Then an optimal penalty value is
chosen by looking for a stabilization in the number of
segments according to the penalty values.
}
\author{
Quentin Grimonprez
}

