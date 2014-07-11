/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2007  Serge Iovleff

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
 * Project:  stkpp::Algebra
 * Purpose:  Define 2D Linear Algebra methods with Real Containers
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_LinAlgebra2D.h
 *  @brief Define Linear Algebra methods for two dimensional
 *  Containers containing Real.
 **/
 
#ifndef STK_LINALGEBRA2D_H
#define STK_LINALGEBRA2D_H

// Point and Vector classes
#include "Arrays/include/STK_Array2DSquare.h"

namespace STK
{
/** @ingroup Algebra
 *  @brief transpose a matrix
 * 
 *  Transpose the Matrix A and give the result in At.
 * 
 *  @param[in] A the matrix to transpose
 *  @param[out] At the transposed matrix
 **/
template < class TContainer2D1, class TContainer2D2>
void transpose( TContainer2D1 const& A, TContainer2D2& At)
{
  // Resize At.
  At.resize(A.cols(), A.rows());
  // copy each column of A in each row of At
  for (int j=A.firstIdxRows(); j<=A.lastIdxRows(); j++)
    for (int i=A.firstIdxCols(); i<=A.lastIdxCols(); i++)
      At(i, j) = A(j, i);
}

/** @ingroup Algebra
 *  @brief Compute the trace of a  square Matrix
 *
 *  @param[in] A the square matrix
 *  @return the trace of the matrix
 **/
template < class TContainer2D >
typename TContainer2D::Type trace( ITContainer<TContainer2D> const& A)
{
  typename TContainer2D::Type sum = 0.0;
  for (int k = A.begin(); k<= A.lastIdx(); k++) sum += A(k, k);
  return sum;
}

} // namespace STK

#endif
// STK_LINALGEBRA2D_H
