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
 * Purpose: implement the AIC criterion.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_AICCriterion.h
 *  @brief In this file we implement the AICCriterion class.
 **/

#include "../../STKernel/include/STK_Exceptions.h"

#include "../include/STK_AICCriterion.h"
namespace STK
{
//Constructor
AICCriterion::AICCriterion( IModelBase const& model)
                          : ICriterion(model)
{}



//Destructor
AICCriterion::~AICCriterion(){}

//* Compute AIC Criterion */
bool AICCriterion::run()
{
  try
  {
    Real loglikelihood    = p_model_->lnLikelihood();
    int freeParameter = p_model_->nbFreeParameters();
    // AIC criteria
    value_  = 2.*(-loglikelihood+freeParameter);
  }
  catch(const Exception& e)
  {
    msg_error_ = e.error();
    return false;
  }
  return true;
}


} // namespace STK



