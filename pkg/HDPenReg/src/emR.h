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
RcppExport SEXP EMlasso(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP lambda, SEXP intercept, SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG);
RcppExport SEXP EMfusedLasso(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP lambda1, SEXP lambda2, SEXP intercept, SEXP maxStep, SEXP burn, SEXP eps, SEXP eps0, SEXP epsCG);
RcppExport SEXP EMlogisticlasso(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP lambda, SEXP intercept, SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG);
RcppExport SEXP EMlogisticfusedLasso(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP lambda1, SEXP lambda2, SEXP intercept, SEXP maxStep, SEXP burn, SEXP eps, SEXP eps0, SEXP epsCG);
RcppExport SEXP EMCVLasso(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP lambda, SEXP nbFolds, SEXP intercept, SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG);
RcppExport SEXP EMCVFusedLasso1D(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP lambda1, SEXP lambda2, SEXP optimL1, SEXP nbFolds, SEXP intercept, SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG);
RcppExport SEXP EMCVFusedLasso2D(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP lambda1, SEXP lambda2, SEXP nbFolds, SEXP intercept, SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG);
RcppExport SEXP EMCVLogisticLasso(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP lambda, SEXP nbFolds, SEXP intercept, SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG);
RcppExport SEXP EMCVLogisticFusedLasso1D(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP lambda1, SEXP lambda2, SEXP optimL1, SEXP nbFolds, SEXP intercept, SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG);
RcppExport SEXP EMCVLogisticFusedLasso2D(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP lambda1, SEXP lambda2, SEXP nbFolds, SEXP intercept, SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG);

#endif
