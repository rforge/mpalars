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

    Contact : Serge.Iovleff@stkpp.org
*/

/*
 * Project:  stkpp::
 * created on: 2 sept. 2013
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_Clust_Util.h
 *  @brief In this file .
 **/


#ifndef STK_CLUST_UTIL_H_
#define STK_CLUST_UTIL_H_

#include "../../STKernel/include/STK_Real.h"

namespace STK
{

namespace Clust
{

/** @ingroup Clustering
 * Default number of iterations in the short runs (used in strategy) */
const int maxIterShortRun = 200;
/**  @ingroup Clustering
 * Default number of iterations in the long run (used in strategy) */
const int maxIterLongRun = 1000;

/**  @ingroup Clustering
 * Default epsilon in the short runs (used in strategy) */
const Real epsilonShortRun = 1e-04;
/**  @ingroup Clustering
 * Default epsilon in the long run (used in strategy) */
const Real epsilonLongRun = Arithmetic<Real>::epsilon();


}  // namespace Clust


}  // namespace STK

#endif /* STK_CLUST_UTIL_H_ */
