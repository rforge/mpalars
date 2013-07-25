/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2013  Serge Iovleff

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  stkpp::Rtools
 * Purpose:  Define tools for converting R data to C data.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/


#ifndef CONVERSION_H
#define CONVERSION_H

#include <Rcpp.h>

namespace Conversion
{
  /** Convert a C matrix in a Rcpp Matrix
   *  @param nbRow number of rows of the matrix
   *  @param nbCol number of column of the matrix
   *  @param matrix the matrix to convert
   *  @return the Rcpp Matrix
   */
  Rcpp::NumericMatrix CMatrixToRcppMatrix(int nbRow, int nbCol, double** matrix);
  
  /** Convert a c vector in a Rcpp Vector
   *  @param dim the dimension of the vector
   *  @param vector the C vector
   *  @return the Rcpp Vector
   */
  Rcpp::NumericMatrix CMatrixToRcppMatrixForInt(int nbRow, int nbCol, int** matrix);

  /** Convert a c vector in a Rcpp Vector
   *  @param dim the dimension of the vector
   *  @param vector the C vector
   *  @return the Rcpp Vector
   */
  Rcpp::NumericVector CVectorToRcppVector(int dim, double* vector);

  /** Convert a c vector in a Rcpp Vector
   *  @param nbCluster number of cluster
   *  @param labels array with the label of the individuals
   *  @return A Rcpp matrix with the partition of in a binary form
   */
  Rcpp::NumericVector CVectorToRcppVectorForInt(int dim, int* vector);
}

#endif /* CONVERSION_H */
