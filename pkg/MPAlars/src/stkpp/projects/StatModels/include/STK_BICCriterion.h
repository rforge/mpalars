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

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  stkpp::Model
 * created on: 22 juil. 2011
 * Purpose: define the BIC criterion.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_BICCriterion.h
 *  @brief In this file we define the BICCriterion class.
 **/

#ifndef STK_BICCRITERION_H
#define STK_BICCRITERION_H

#include "STK_ICriterion.h"

namespace STK
{
/** @ingroup StatModels
 *  @brief Derived class of Criterion for computing the BIC Criterion
 *  The Bic criteria is a penalization of the likelihood given by the formula
 *  \f[
 *  -2 \cdot \ln{p(x|k)} \approx \mathrm{BIC} = {-2 \cdot \ln{L} + D \ln(n) }
 *  \f]
 *  where \f$ L \f$ represents the likelihood of the observations, \f$ D \f$ the
 *  number of free parameter of the model and \f$ n \f$ the number of sample.
 **/
class BICCriterion : public ICriterion
{
  public:
    /** Constructor.
     *  @param model the model to evaluate the criterion
     **/
    BICCriterion( IModelBase const& model);
    /** virtual destructor. */
    virtual ~BICCriterion();
    /** implementation of the virtual method run */
    virtual bool run();
};

} // namespace STK

#endif
