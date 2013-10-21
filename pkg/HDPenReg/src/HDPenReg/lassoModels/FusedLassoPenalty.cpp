/*--------------------------------------------------------------------*/
/*     Copyright (C) 2013-2013  Serge Iovleff, Quentin Grimonprez

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

    Contact : quentin.grimonprez@inria.fr
*/

/*
 * Project:  MPAGenomics::
 * created on: 13 sept. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file FusedLassoPenalty.cpp
 *  @brief In this file, implementation of the method of the @c FusedLassoPenalty class .
 **/

#include "FusedLassoPenalty.h"

namespace HD
{

  /*default constructor*/
  FusedLassoPenalty::FusedLassoPenalty()
  : lambda1_(), lambda2_(), mainDiagonal_(), offDiagonal_(), sigma2_(1.), eps_(STK::Arithmetic<STK::Real>::epsilon())
  {
  }

  /* Constructor
   *  @param lambda1 penalization parameter for the l1-norm of the estimates
   *  @param lambda2 penalization parameter for the l1-norm of the difference between successive estimates
   *  @param eps epsilon to add to denominator of fraction to avoid zeros.
   */
  FusedLassoPenalty::FusedLassoPenalty(STK::Real lambda1, STK::Real lambda2, STK::Real eps)
                                      : lambda1_(lambda1)
                                      , lambda2_(lambda2)
                                      , mainDiagonal_()
                                      , offDiagonal_()
                                      , sigma2_(1)
                                      , eps_(eps)
  {
  }


  /*
   * Copy constructor
   * @param penalty LassoPenalty object to copy
   */
  FusedLassoPenalty::FusedLassoPenalty(FusedLassoPenalty const& penalty)
              : IPenalty(penalty)
              , lambda1_(penalty.lambda1())
              , lambda2_(penalty.lambda2())
              , mainDiagonal_(penalty.mainDiagonal())
              , offDiagonal_(penalty.offDiagonal())
              , sigma2_(penalty.sigma2())
              , eps_(penalty.eps())
  {
  }

  /*
   * clone
   */
  FusedLassoPenalty* FusedLassoPenalty::FusedLassoPenalty::clone() const
  {
    return new FusedLassoPenalty(*this);
  }

  /*
   * @param beta current estimates
   * @return t(beta) * matrixB * beta
   */
  STK::Real FusedLassoPenalty::penaltyTerm(STK::CVectorX const& beta) const
  {
    STK::Real pen(0);
//    pen = beta.dot(matrixB_ * beta);
    //check the size
    if(beta.sizeRows() != mainDiagonal_.sizeRows())
      throw(STK::out_of_range("size mismatch."));

    if(beta.sizeRows() > 0)
    {
      if(beta.sizeRows() == 1)
        pen = beta[1] * mainDiagonal_[1] * beta[1];
      else
      {
        pen = beta[1] * (mainDiagonal_[1] * beta[1] +  offDiagonal_[1] * beta[2]);
        if(beta.sizeRows() > 2)
        {
          for(int i = 2; i < beta.sizeRows(); i++)
            pen += beta[i] * ( offDiagonal_[i-1] * beta[i-1] + mainDiagonal_[i] * beta[i] +  offDiagonal_[i] * beta[i+1] );
        }
        pen += beta.back() * (beta.back() * mainDiagonal_.back() + beta[beta.sizeRows()-1] * offDiagonal_[beta.sizeRows()-1] );
      }

    }

    return pen;
  }

  /*
   * update sigma2 and the fused lasso penalty
   * @param beta current estimates
   * @param normResidual ||y-X*beta||_2^2
   */
  void FusedLassoPenalty::update(STK::CVectorX const& beta)
  {
    updatePenalty(beta);
  }

  /* update the penalty matrix : matrixB_
   *  @param beta current estimates
   */
  void FusedLassoPenalty::updatePenalty(STK::CVectorX const& beta)
  {
    //resize to the current size
    offDiagonal_.resize(beta.sizeRows()-1);
    mainDiagonal_.resize(beta.sizeRows());

    if(beta.sizeRows() == 1)
      mainDiagonal_[1] = lambda1_/(std::abs(beta[1]) + eps_);
    else
    {
      offDiagonal_[1] = -lambda2_/(std::abs(beta[2]-beta[1]) + eps_);
      mainDiagonal_[1] = lambda1_/(std::abs(beta[1]) + eps_) - offDiagonal_[1];

      if(beta.sizeRows() > 2)
      {
        for(int i = 2; i < beta.sizeRows(); i++)
        {
          offDiagonal_[i] = -lambda2_/(std::abs(beta[i+1]-beta[i]) + eps_);
          mainDiagonal_[i] = lambda1_/(std::abs(beta[i]) + eps_) - offDiagonal_[i] - offDiagonal_[i-1];
        }
      }

      mainDiagonal_[beta.sizeRows()] = 1/(std::abs(beta[beta.sizeRows()]) + eps_) - offDiagonal_[beta.sizeRows()-1];
    }
  }

}


