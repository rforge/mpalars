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
 * Purpose: implement the BIC criterion.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_BICCriterion.h
 *  @brief In this file we implement the BICCriterion class.
 **/
#include "../../STKernel/include/STK_Exceptions.h"

#include "../include/STK_BICCriterion.h"

namespace STK
{
//Constructor
BICCriterion::BICCriterion( IModelBase const& model)
                          : ICriterion(model)
{}



//Destructor
BICCriterion::~BICCriterion(){}

//* Compute BIC Criterion */
bool BICCriterion::run()
{
  try
  {
    // BIC criteria
    value_  = (-2.*p_model_->lnLikelihood())+(p_model_->nbFreeParameter()*p_model_->lnNbSample());
  }
  catch( Exception const& e)
  {
    msg_error_ = e.error();
    return false;
  }
  return true;
}


} // namespace STK



