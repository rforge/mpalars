#ifndef _EMR_H
#define _EMR_H

#include <Rcpp.h>
#include <vector>
#include "stkpp/projects/STKernel/include/STK_Real.h"
#include "stkpp/include/STKernel.h"
#include "stkpp/include/Arrays.h"

/*
 * note : RcppExport is an alias to `extern "C"` defined by Rcpp.
 *
 * It gives C calling convention to the rcpp_hello_world function so that 
 * it can be called from .Call in R. Otherwise, the C++ compiler mangles the 
 * name of the function and .Call can't find it.
 *
 * It is only useful to use RcppExport when the function is intended to be called
 * by .Call. See the thread http://thread.gmane.org/gmane.comp.lang.r.rcpp/649/focus=672
 * on Rcpp-devel for a misuse of RcppExport
 */


void copySTKVectorInSTDVector(STK::CVectorX const& stkvector, std::vector<STK::Real> &stdvector);
void copySTKArray2DVectorInSTDVector(STK::Array2DVector<int> const& stkvector, std::vector<int> &stdvector);
RcppExport SEXP EMlasso(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP lambda, SEXP maxStep, SEXP eps);


#endif
