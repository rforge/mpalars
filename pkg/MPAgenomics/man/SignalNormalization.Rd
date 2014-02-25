\name{SignalNormalization}
\alias{SignalNormalization}
\title{normalization process.}
\usage{
SignalNormalization(dataFolder, chipType, normalTumorArray,
  genotypeCallsMethod = "naive", savePlot = TRUE, tags = NULL)
}
\arguments{
  \item{dataFolder}{name of the data set.}

  \item{chipType}{type of the chip used for the data.}

  \item{normalTumorArray}{only if you have normal and tumor
  profile in your data folder. A csv file or a data.frame
  with 2 columns: "normal" and "tumor". The first column
  contains the name of normal files and the second the
  names of associated tumor files.}

  \item{genotypeCallsMethod}{method used for genotypage,
  default is "naive".}

  \item{savePlot}{if TRUE, graphics of the CN signal and
  allele B fraction signal will be saved in the figures
  folder.}

  \item{tags}{appear in the different file names (cdf, ugp,
  ufl) of the chip. For no tag, use tags=NULL (default =
  NULL). See details for more information.}
}
\description{
normalization process for estimating raw copy-numbers and
allele B fraction.
}
\details{
You have to respect the aroma architecture: <working
directory> +- annotationData/ | +- chipTypes/ | +-
<chipType>/ <-- must match exactly the name of the CDF file
(fullname minus tags) | +- CDF file(s) and other annotation
(possibly subdirectories) | +- rawData/ +- <nameOfDataSet>/
+- <chipType>/ <-- must match exactly a chip type folder
under annotationData/ +- CEL files

All the cdf chip file names must follow the following rule
: <chipType>,<Tags>.cdf Multiples tags must be separated by
a comma. If there is no tags, the pattern is <chipType>.cdf
}
\author{
Quentin Grimonprez
}

