#' @useDynLib HDPenReg
#' @import Rcpp
#' @import methods
#' 
#' @title Algorithms for lasso and fused-lasso problems.
#' @docType package
#' @aliases HDPenReg-package, HDPenReg
#' @name HDPenReg-package
#' @description This package contains algorithms for lasso and fused-lasso problems.
#' It contains an implementation of the lars algorithm [1],
#' for the lasso and fusion penalization and EM-based algorithms for (logistic) lasso and fused-lasso.
#' 
#' @details
#' 
#'   \tabular{ll}{
#' Package: \tab HDPenReg\cr
#' Type: \tab Package\cr
#' Version: \tab 0.89.5\cr
#' Date: \tab 2014-06-24\cr
#' License: \tab GPL (>=2) \cr
#' }
#' 
#' The main function is \link{HDlars}. 
#' 
#' 
#' @author Maintainer: Quentin Grimonprez <quentin.grimonprez@@inria.fr>
#' 
#' @examples 
#' #see vignette
#' #vignette("HDPenReg")
#' 
#'  
#' @keywords package
NULL