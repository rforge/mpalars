#ifndef _LARSR_H
#define _LARSR_H

#include <Rcpp.h>
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

void convertToVector(SEXP const& rVector, SEXP size, STK::CVectorX &output);
void convertToArray(SEXP const& rMatrix, STK::CArrayXX &output);
RcppExport SEXP lars(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP maxStep, SEXP intercept, SEXP eps);
RcppExport SEXP cvlars(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP maxStep, SEXP intercept, SEXP eps, SEXP nbFold, SEXP partition, SEXP index);
RcppExport SEXP fusion(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP maxStep, SEXP intercept, SEXP eps);

#endif
