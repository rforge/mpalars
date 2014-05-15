#' @import R.utils snowfall changepoint glmnet cghseg
#' 
#' @title Multi Patients Analysis of Genomic Markers
#' @docType package
#' @aliases MPAgenomics-package, MPAgenomics
#' @name MPAgenomics-package
#' @description  This package provides functions to preprocess and analyze genomic data. 
#' The package was initially developped to select genomic markers associated with a given phenotype when several samples are available.
#' In this context, markers refer to SNPs or copy number variations which are designed on the arrays. 
#'
#' The package also enables to preprocess all samples individually in order to keep maximum information from the original signals and improve the multi-patients analysis. In particular, this is useful to keep quantitative data for SNPs rather than usual genotype calls (AA, AB or BB) when these states are not relevant (eg in cancer studies where the number of copies differs from two copies).
#'
#' 
#' @details
#' 
#'   \tabular{ll}{
#' Package: \tab MPAgenomics\cr
#' Type: \tab Package\cr
#' Version: \tab 1.0.3\cr
#' Date: \tab 2014-05-15\cr
#' License: \tab GPL (>=2) \cr
#' }
#' 
#' 
#' 
#' @author Quentin Grimonprez with contributions from Guillemette Marot and Samuel Blanck
#'
#' Maintainer: Samuel Blanck <samuel.blanck@@inria.fr>
#'  
#' 
#' @examples 
#' see the vignette for detailed examples
#' vignette("MPAgenomics")
#'  
#' @keywords package
NULL
