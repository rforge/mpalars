#include "larsR.h"
#include "HDPenReg/lars/Lars.h"
#include "HDPenReg/lars/Cvlars.h"
#include "HDPenReg/lars/Fusion.h"

#include <iostream>

using namespace Rcpp;

RcppExport SEXP lars(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP maxStep, SEXP intercept, SEXP eps)
{ return larsmain(data, response, nbIndiv, nbVar, maxStep, intercept, eps);}

RcppExport SEXP fusion(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP maxStep, SEXP intercept, SEXP eps)
{ return fusionmain(data, response, nbIndiv, nbVar, maxStep, intercept, eps);}

RcppExport SEXP cvlars(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP maxStep, SEXP intercept, SEXP eps, SEXP nbFold, SEXP partition, SEXP index, SEXP mode)
{ return cvlarsmain(data, response, nbIndiv, nbVar, maxStep, intercept, eps, nbFold, partition, index, mode);}
