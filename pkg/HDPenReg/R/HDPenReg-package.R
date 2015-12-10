#' @useDynLib HDPenReg
#' @import rtkore
#' @import methods
#' @importFrom graphics abline axis lines points
#' @importFrom stats rbeta rbinom rpois
#' 
#' @title Algorithms for lasso and fused-lasso problems.
#' @docType package
#' @aliases HDPenReg-package HDPenReg
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
#' Version: \tab 0.93\cr
#' Date: \tab 2015-10-09\cr
#' License: \tab GPL (>=2) \cr
#' }
#' 
#' The main function is \code{\link{HDlars}}, it performs the lars algorithm [1] for solving lasso problem. Some functions are provided to perform
#'  cross-validation and prediction. An implementation of the lasso and fused-lasso with an EM algorithm is also provided.
#' 
#' 
#' @source [1] Efron, Hastie, Johnstone and Tibshirani (2003) "Least Angle Regression" (with discussion) Annals of Statistics
#' 
#' @author Maintainer: Quentin Grimonprez <quentin.grimonprez@@inria.fr>
#' 
#' @examples 
#' \dontrun{
#' #see vignette
#' vignette("HDPenReg")
#' }
#'  
#'  @seealso \code{\link{HDlars}} \code{\link{HDcvlars}}
#'  
#' @keywords package
NULL
