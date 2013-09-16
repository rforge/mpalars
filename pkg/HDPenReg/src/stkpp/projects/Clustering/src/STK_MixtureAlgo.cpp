/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2012  Serge Iovleff

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
 * Project:  stkpp::Clustering
 * created on: 16 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 * Originally created by Parmeet bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_IMixtureModel.cpp
 *  @brief In this file we implement the run method of the mixture algorithms.
 **/

#include "../include/STK_MixtureAlgo.h"
#include "../include/STK_IMixtureModelBase.h"

namespace STK
{

bool CEMAlgo::run()
{
  try
  {
    Real currentLikelihood = -STK::Arithmetic<Real>::max();
    for (int iter = 0; iter < nbIterMax_; ++iter)
    {
      p_model_->ceStep();
      p_model_->mStep();
      p_model_->computeLnLikelihood();
      // no abs as the likelihood should increase
      if ( (p_model_->lnLikelihood() - currentLikelihood) < epsilon_) break;
      currentLikelihood = p_model_->lnLikelihood();
    }
  }
  catch (Exception const& e)
  {
    msg_error_ = e.error();
    return false;
  }
  return true;
}

bool EMAlgo::run()
{
  try
  {
    for (int iter = 0; iter < this->nbIterMax_; ++iter)
    {
      Real currentLikelihood = -STK::Arithmetic<Real>::max();
      for (int iter = 0; iter < nbIterMax_; ++iter)
      {
        p_model_->eStep();
        p_model_->mStep();
        p_model_->computeLnLikelihood();
        // no abs as the likelihood should increase
        if ( (p_model_->lnLikelihood() - currentLikelihood) < epsilon_) break;
        currentLikelihood = p_model_->lnLikelihood();
      }
      // compute zi
      p_model_->mapStep();
    }
  }
  catch (Exception const& e)
  {
    msg_error_ = e.error();
    return false;
  }
  return true;
}

bool SEMAlgo::run()
{
  try
  {
    for (int iter = 0; iter < this->nbIterMax_; ++iter)
    {
      p_model_->seStep();
      p_model_->mStep();
    }
    p_model_->computeLnLikelihood();
  }
  catch (Exception const& e)
  {
    msg_error_ = e.error();
    return false;
  }
  return true;
}

} // namespace STK
