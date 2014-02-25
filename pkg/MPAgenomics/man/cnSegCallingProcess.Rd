\name{cnSegCallingProcess}
\alias{cnSegCallingProcess}
\title{Segment a copy-number signal and call the found segments.}
\usage{
cnSegCallingProcess(dataSetName, normalTumorArray, chromosome = 1:22,
  Lambda = NULL, listOfFiles = NULL, onlySNP = TRUE, savePlot = TRUE,
  nclass = 3, cellularity = 1, ...)
}
\arguments{
  \item{dataSetName}{name of the data-set folder in the
  rawData folder containing the signals to used.}

  \item{normalTumorArray}{only if you have normal and tumor
  profile in your data folder. A csv file or a data.frame
  with 2 columns: "normal" and "tumor". The first column
  contains the name of normal files and the second the
  names of associated tumor files.}

  \item{chromosome}{A vector containing the chromosomes to
  segment.}

  \item{Lambda}{vector containing all the penalization
  values to test for the segmentation. If no values are
  provided, default values will be used.}

  \item{listOfFiles}{A vector containing the names of the
  files from the dataSetName to use.}

  \item{onlySNP}{If TRUE, only the SNP probes will be
  used.}

  \item{savePlot}{If TRUE, save the segmented signal in
  figures folder.}

  \item{nclass}{The number of levels to be used for
  calling. Either 3 (loss, normal, gain), 4 (including
  amplifications), 5 (including double deletions).}

  \item{cellularity}{Percentage of tumored cells in the
  sample}

  \item{...}{Other parameters of CGHcall function}
}
\value{
a data.frame containg columns : \describe{
\item{sampleNames}{Name of the file.} \item{chrom}{The
chromosome of the segment.} \item{chromStart}{The starting
position (in bp) of a segment. This position is not
included in the segment.} \item{chromEnd}{The ending
position (in bp) of a segment.This position is included in
the segment.} \item{probes}{Number of probes in the
segment.} \item{means}{Mean of the segment.}
\item{calls}{The calling of segment ("double loss", "loss",
"normal", "gain" or "amplification").} }
}
\description{
This function applies the PELT method for segment each
signal of the dataset and launches CGHcall for calling
segments and detect aberrations. Results will be stored in
a BED file in the segmentation folder of the aroma
architecture.
}
\examples{
#DO NOT EXECUTE
# seg1=cnSegCallingProcess("data1",normalTumorArray,chromosome=20:21)
# seg2=cnSegCallingProcess("data2",chromosome=20:21)
}
\author{
Quentin Grimonprez
}

