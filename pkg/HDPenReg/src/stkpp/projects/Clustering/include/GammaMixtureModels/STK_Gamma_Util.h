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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project: stkpp::Clustering
 * created on: 5 sept. 2013
 * Author:  iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_Gamma_Util.h
 *  @brief In this file we implement some functionalities used in all gamma models.
 **/

#ifndef GAMMA_UTIL_H
#define GAMMA_UTIL_H

#include "../../../Arrays/include/STK_Array2D.h"
#include "STK_GammaParameters.h"
#include "STK_GammaComponent.h"

namespace STK
{
/** @ingroup Clustering
 *  enclosing class for the methods common to all the gamma models
 **/
template<class Component>
struct GammaUtil
{
  typedef typename Component::Array Array;
  typedef typename Component::Parameters Parameters;
  typedef typename Array::Col ColVector;

  /** compute the moment.
   *  @param components the components with the parameters to initialize
   *  @param p_tik the tik
   **/
  static void moments(Array1D< Component* >& components, Array2D<Real> const* p_tik)
  {
    if (components.size() <= 0) return;
    for (int k= p_tik->firstIdxCols(); k <= p_tik->lastIdxCols(); ++k)
    {
      Parameters* paramk = components[k]->p_param();
      Array const* p_data = components[k]->p_data();
      ColVector tik(p_tik->col(k), true); // create a reference

      for (int j=p_data->firstIdxCols(); j<=p_data->lastIdxCols(); ++j)
      {
        Real mean =  p_data->col(j).wmeanSafe(tik);
        if ((mean<=0)||Arithmetic<Real>::isNA(mean)) throw Clust::estimFail_;
        paramk->mean_[j] = mean;
        Real meanLog =  p_data->col(j).log().wmeanSafe(tik);
        if ((meanLog<=0)||Arithmetic<Real>::isNA(meanLog)) throw Clust::estimFail_;
        paramk->meanLog_[j] = meanLog;
        Real variance =  p_data->col(j).wvarianceSafe(tik);
        if ((variance<=0)||Arithmetic<Real>::isNA(variance)) throw Clust::estimFail_;
        paramk->variance_[j] = variance;
      }
    }
  }
};

}  // namespace STK

#endif /* GAMMA_UTIL_H */
