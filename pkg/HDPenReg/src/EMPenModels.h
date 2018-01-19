#ifndef EMPENMODELS_H
#define EMPENMODELS_H

#include <RTKpp.h>

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

RcppExport SEXP EMlassoMain( SEXP data, SEXP response
                           , SEXP lambda, SEXP intercept
                           , SEXP maxStep, SEXP burn
                           , SEXP threshold, SEXP eps, SEXP epsCG);

RcppExport SEXP EMlogisticLassoMain( SEXP data, SEXP response
                                   , SEXP lambda, SEXP intercept
                                   , SEXP maxStep, SEXP burn
                                   , SEXP threshold, SEXP eps, SEXP epsCG);

RcppExport SEXP EMfusedLassoMain( SEXP data, SEXP response
                                , SEXP lambda1, SEXP lambda2, SEXP intercept
                                , SEXP maxStep, SEXP burn
                                , SEXP eps, SEXP eps0, SEXP epsCG);

RcppExport SEXP EMlogisticFusedLassoMain( SEXP data, SEXP response
                                        , SEXP lambda1, SEXP lambda2, SEXP intercept
                                        , SEXP maxStep, SEXP burn
                                        , SEXP eps, SEXP eps0, SEXP epsCG);

RcppExport SEXP cvEMlassoMain( SEXP data, SEXP response
                             , SEXP lambda
                             , SEXP nbFolds
                             , SEXP intercept
                             , SEXP maxStep, SEXP burn
                             , SEXP threshold, SEXP eps, SEXP epsCG);

RcppExport SEXP cvEMfusedLasso1DMain( SEXP data, SEXP response
                                    , SEXP lambda1, SEXP lambda2
                                    , SEXP optimL1, SEXP nbFolds
                                    , SEXP intercept, SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG);

RcppExport SEXP cvEMfusedLasso2DMain( SEXP data, SEXP response
                                    , SEXP lambda1, SEXP lambda2
                                    , SEXP nbFolds
                                    , SEXP intercept, SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG);

RcppExport SEXP cvEMlogisticLassoMain( SEXP data, SEXP response
                                     , SEXP lambda, SEXP nbFolds
                                     , SEXP intercept, SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG);

RcppExport SEXP cvEMlogisticFusedLasso1DMain( SEXP data, SEXP response
                                            , SEXP lambda1, SEXP lambda2
                                            , SEXP optimL1, SEXP nbFolds
                                            , SEXP intercept, SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG);

RcppExport SEXP cvEMlogisticFusedLasso2DMain( SEXP data, SEXP response
                                                , SEXP lambda1, SEXP lambda2
                                                , SEXP nbFolds
                                                , SEXP intercept
                                                , SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG);



#endif
