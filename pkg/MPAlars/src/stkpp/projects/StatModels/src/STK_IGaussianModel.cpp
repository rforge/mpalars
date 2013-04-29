/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2011  Serge Iovleff

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

 Contact : Serge.Iovleff@stkpp.org
 */

/*
 * Project:  stkpp::
 * created on: 13 août 2011
 * Purpose: implemen the interface class IGaussianModel .
 * Author:   iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_IGaussianModel.cpp
 *  @brief In this file we implement the Interface class IGaussianModel.
 **/

#include<cmath>

#ifdef STK_DEBUG
#include "../../../include/DManager.h"
#endif

#include "../include/STK_IGaussianModel.h"

#include "../../STatistiK/include/STK_Stat_UnivariateReal.h"
#include "../../STKernel/include/STK_Arithmetic.h"

namespace STK
{

/**
 * Compute the gaussian log likehood of a one dimensionnal gaussian model.
 * @param data the data set
 * @param mu the mean of the gaussian law
 * @param sigma the variance of the gaussian law
 * @return
 */
Real univariateGaussianLnLikelihood(Vector const& data, Real const& mu, Real const& sigma)
{
  int first = data.firstIdx(), last = data.lastIdx(), nbSample = data.size();
  if (sigma)
  {
    Real scale = 0., std = std::sqrt((double)sigma);
    // compute scale
    for (int i= first; i <= last; ++i)
    {
      scale = std::max(scale, std::abs((data[i]-mu)/std));
    }
    Real sum = 0;
    if (scale)
    {
      // compute sum of centered variable
      for (int i= first; i <= last; ++i)
      {
        Real res = ((data[i]-mu)/std)/scale;
        sum += (res*res);
      }
    }
   return - (0.5*sum*scale*scale + nbSample * (std::log((double)std) + Const::_LNSQRT2PI_));
  }
  // 0 variance
  return -STK::Arithmetic<Real>::infinity();
}
/**
 * Compute the gaussian log likelihood of a diagonal gaussian model.
 * @param data the data set
 * @param mu the mean of the gaussian law
 * @param sigma the (diagonal) covairance matrix
 * @return
 */
Real diagonalGaussianLnLikelihood(Matrix const& data, Point const& mu, MatrixSquare const& sigma)
{
  int first = mu.firstIdx(), last = mu.lastIdx();
  // compute for each row tjhe gaussian ln-likehood
  Real sum = 0.;
  for (int j = first; j<= last; ++j)
  {
    //stk_cout << data.col(j) << _T("\n");
    sum +=univariateGaussianLnLikelihood(data.col(j), mu[j], sigma(j,j));
  }
  return sum;
}

} // namespace STK
