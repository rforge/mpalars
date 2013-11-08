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
 * created on: 21 juin 2013
 * Author:   Quentin Grimonprez
 **/

/** @file EnetPenalty.cpp
 *  @brief In this file, methods of class @c EnetPenalty are implemented.
 **/



#include "EnetPenalty.h"

namespace HD
{
  /* Constructor
   *  @param lambda1 penalization parameter for the l1-norm of the estimates
   *  @param lambda2 penalization parameter for the l2-norm of the estimates
   *  @param n size of sample
   *  @param p size of Penalty (number of covariates)
   */
  EnetPenalty::EnetPenalty( STK::Real lambda1, STK::Real lambda2, int n, int p)
                          : IPenalty()
                          , lambda1_(lambda1)
                          , lambda2_(lambda2)
                          , sqrtInvPenalty_(p,0)
                          , sigma2_(1)
                          , n_(n)
                          , p_(p)
  {
  }


  /*
   * Copy constructor
   * @param penalty LassoPenalty object to copy
   */
  EnetPenalty::EnetPenalty(EnetPenalty const& penalty)
              : IPenalty(penalty)
              , lambda1_(penalty.lamba1())
              , lambda2_(penalty.lamba2())
              , sqrtInvPenalty_(penalty.sqrtInvPenalty())
              , sigma2_(penalty.sigma2())
              , n_(penalty.n())
              , p_(penalty.p())
  {
  }

  /*clone*/
  //EnetPenalty* EnetPenalty::clone() const
  //{
  //  return new EnetPenalty(*this);
  //}

  /*
   * update sigma2 and the lasso penalty
   * @param beta current estimates
   * @param normResidual ||y-X*beta||_2^2
   */
  void EnetPenalty::update(STK::CVectorX const& beta, STK::Real const& normResidual)
  {
    updatePenalty(beta);
    updateSigma2(beta,normResidual);
  }

  void EnetPenalty::update(STK::CVectorX const& beta)
  {
    updatePenalty(beta);
  }


  /*
   * @param x a vector of length p_
   * @return the product invPenalty_*x
   */
  STK::CVectorX EnetPenalty::multInvPenalty(STK::CVectorX const& x) const
  {
    STK::CVectorX a;
    a = sqrtInvPenalty_.square() * x;

    return a;
  }

  /*
   * @param x a vector of length p_
   * @return the product invPenalty_.sqrt()*x
   */
  STK::CVectorX EnetPenalty::multSqrtInvPenalty(STK::CVectorX const& x) const
  {
    STK::CVectorX a;
    a = sqrtInvPenalty_ * x;

    return a;
  }

  /* penalty term
   *  @param beta current estimates
   *  @return t(beta) * penalty * beta
   */
  STK::Real EnetPenalty::penaltyTerm(STK::CVectorX const& beta) const
  {
    STK::Real penaltyTerm = beta.dot( (sqrtInvPenalty_.square()).inverse() * beta);


    return penaltyTerm;
  }


  /*
   * update the penalty
   * @param beta current estimates
   */
  void EnetPenalty::updatePenalty(STK::CVectorX const& beta)
  {
    sqrtInvPenalty_.resize(beta.sizeRows());

    for(int i = 1; i <= beta.sizeRows(); i++)
      sqrtInvPenalty_[i] = std::sqrt(std::abs(beta[i]) / ( lambda1_ + lambda2_ * std::abs(beta[i]) ) );
  }

  /*
   * update sigma2
   * @param beta current estimates
   * @param normResidual ||y-X*beta||_2^2
   */
  void EnetPenalty::updateSigma2(STK::CVectorX const& beta, STK::Real const& normResidual)
  {
    STK::Real normBeta=beta.abs().sum() * lambda1_;
    sigma2_ = (normResidual + normBeta + beta.dot(beta) ) / (n_ + p_ - 1.);
  }
}
