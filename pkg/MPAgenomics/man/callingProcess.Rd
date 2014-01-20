\name{callingProcess}
\alias{callingProcess}
\title{Calling aberrations in segmented copy-number signal.}
\usage{
  callingProcess(segmentData, nclass = 5, cellularity = 1,
    verbose = TRUE, ...)
}
\arguments{
  \item{segmentData}{A list (see details).}

  \item{nclass}{The number of levels to be used for
  calling. Either 3 (loss, normal, gain), 4 (including
  amplifications), 5 (including double deletions).}

  \item{cellularity}{A vector of cellularities ranging from
  0 to 1 to define the contamination of your sample with
  healthy cells (1 = no contamination).}

  \item{verbose}{If TRUE, print some details.}

  \item{...}{other options of CGHcall functions}
}
\value{
  A list with the same element as segmentData list and
  \describe{ \item{calls}{A matrix, of the same size as
  inputSegmented$copynumber matrix, containing the label of
  each point.  -2=double loss, -1=loss, 0=normal, 1=gain,
  2=amplification.} \item{segment}{A data.frame that sums
  up the different segment found.} \item{probdloss}{(if you
  have ran CGHcall with nclass=5) A matrix of the same size
  as inputSegmented$copynumber matrix. It contains the
  probability for each segmented copynumber to be a double
  loss.} \item{probloss}{A matrix of the same size as
  inputSegmented$copynumber matrix. It contains the
  probability for each segmented copynumber to be a loss.}
  \item{probdnorm}{A matrix of the same size as
  inputSegmented$copynumber matrix. It contains the
  probability for each segmented copynumber to be normal.}
  \item{probdgain}{A matrix of the same size as
  inputSegmented$copynumber matrix. It contains the
  probability for each segmented copynumber to be a gain.}
  \item{probdamp}{(if you have ran CGHcall with nclass=4 or
  5) A matrix of the same size as inputSegmented$copynumber
  matrix. It contains the probability for each segmented
  copynumber to be an amplification.} }
}
\description{
  Launch the process of segmentation labeling. This
  function uses dunctions from CGHcall package developped
  by Sjoerd Vosse, Mark van de Wiel and Ilari Scheinin. See
  the CGHcall package for more details.
}
\details{
  segmentData is a list containing: \describe{
  \item{copynumber}{A matrix. Each column contains a signal
  of copynumber for a profile. Each row corresponds to a
  genomic position of a probe.} \item{segmented}{A matrix
  of the same size as copynumber. It contains the segmented
  signals.} \item{chromosome}{A vector of length
  nrow(copynumber) containing the number of the chromosome
  for each position.} \item{startPos}{A vector of length
  nrow(copynumber) containing the starting genomic position
  of each probe.} \item{featureNames}{A vector of length
  nrow(copynumber) containing the names of each probes.}
  \item{sampleNames}{A vector of length ncol(copynumber)
  containing the names of each profiles.} }
}

