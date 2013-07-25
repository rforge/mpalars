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

#include "Conversion.h"

namespace Conversion
{


//To spin off double** in Rcpp::NumericMatrix
Rcpp::NumericMatrix CMatrixToRcppMatrix(int nbRow, int nbCol, double** matrix)
{
  //NumericMatrix for matrix
  Rcpp::NumericMatrix matrixOutput(nbRow,nbCol);
  for(int i=0;i<nbRow;i++)
  {
    for(int j=0;j<nbCol;j++)
    {
      matrixOutput(i,j) = matrix[i][j];
    }
  }
  return matrixOutput;
}


//To spin off int** in Rcpp::NumericMatrix
Rcpp::NumericMatrix CMatrixToRcppMatrixForInt(int nbRow, int nbCol, int** matrix)
{
    //NumericMatrix for matrix
    Rcpp::NumericMatrix matrixOutput(nbRow,nbCol);
    for(int i=0;i<nbRow;i++)
    {
      for(int j=0;j<nbCol;j++)
      {
        matrixOutput(i,j) = matrix[i][j];
      }
    }
    return matrixOutput;
}

//To spin off double* in Rcpp::NumericVector
Rcpp::NumericVector CVectorToRcppVector(int dim, double* vector)
{
  // NumericVector for vector
  Rcpp::NumericVector vectorOutput(dim);
  for(int i=0;i<dim;i++)
  {
      vectorOutput(i) = vector[i];
  }
  return vectorOutput;
}

//To spin off int* in Rcpp::NumericVector
Rcpp::NumericVector CVectorToRcppVectorForInt(int dim, int* vector)
{
  // NumericVector for vector
  Rcpp::NumericVector vectorOutput(dim);
  for(int i=0;i<dim;i++)
  {
      vectorOutput(i) = vector[i];
  }
  return vectorOutput;
}

}


















