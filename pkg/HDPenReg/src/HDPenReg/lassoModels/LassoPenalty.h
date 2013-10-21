/*--------------------------------------------------------------------*/
/*     Copyright (C) 2013-2013  Serge Iovleff, Quentin Grimonprez

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public
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
 * created on: 28 mai 2013
 * Author:   Quentin Grimonprez
 **/

/** @file LassoPenalty.h
 *  @brief In this file, definition of the class lasso penalty and the functor lassoMultiplicator.
 **/


#ifndef LASSOPENALTY_H_
#define LASSOPENALTY_H_

#include "PenalizedModels.h"
#include "IPenalty.h"

namespace HD
{
  /**functor for CG  */
  struct LassoMultiplicator
  {
      /**
       * functor for the CG. A=sigma2*I+invPenalty.sqrt()*tX*X*invPenalty.sqrt()
       * @param x vector of length sizeCols(A)
       * @return A*x
       */
      STK::CVectorX operator()(STK::CVectorX &x) const
      {
        STK::CVectorX a(x.sizeRows());
        //a = sigI*x+invD*tX*X*invD*x
        a = (*p_sigma2_ * x) + ((*p_invPenalty_ * p_data_->transpose()) * (((*p_data_) * (*p_invPenalty_ * x))));

        return   a ;
      }

      /**
       * Constructor of the functor
       * @param p_data constant pointer on the data
       * @param p_invPenalty constant pointer on the current estimates of invPenalty
       * @param p_sigma2 constant pointer on the current estimates of sigma2
       */
      LassoMultiplicator(STK::CArrayXX const* p_data = 0, STK::Array2DDiagonal<STK::Real> const* p_invPenalty = 0, STK::Real const* p_sigma2 = 0)
                        : p_data_(p_data), p_invPenalty_(p_invPenalty), p_sigma2_(p_sigma2)
      {
      }

      ///pointer to the current data
      STK::CArrayXX const* p_data_;
      ///pointer to the penalty matrix
      STK::Array2DDiagonal<STK::Real> const* p_invPenalty_;
      ///matrix to sigma2
      STK::Real const* p_sigma2_;
  };


  /** @ingroup lassoModels
   *  @brief The class LassoPenalty derived from the @c IPenalty class.
   *  It contains the matrix penalty associated to the lasso problem.
   */
  class LassoPenalty : public IPenalty
  {
    public:
      /**default constructor*/
      LassoPenalty();
      /** Constructor
       *  @param lambda penalization parameter for the l1-norm of the estimates
       *  @param n size of sample
       *  @param p size of Penalty (number of covariates)
       */
      LassoPenalty(STK::Real lambda);

      /** Copy constructor
       *  @param penalty LassoPenalty object to copy
       */
      LassoPenalty(LassoPenalty const& penalty);

      /** destructor */
      inline virtual ~LassoPenalty() {};

      /**clone*/
      LassoPenalty* clone() const;

      //getter
      /**@return lambda parameter of the lasso */
      inline STK::Real const& lambda() const {return lambda_;}
      /**@return invPenalty diagonal matrix containing |beta_i| / lambda */
      inline STK::Array2DDiagonal<STK::Real> const& invPenalty() const {return sqrtInvPenalty_;}
      /**@return sigma2 variance of the response*/
      inline STK::Real const& sigma2() const { return sigma2_;}
      /**@return A constant pointer on the matrix penalty*/
      inline STK::Array2DDiagonal<STK::Real> const* p_invPenalty() const {return &sqrtInvPenalty_;}
      /**@return A constant pointer on sigma2*/
      inline STK::Real const*  p_sigma2() const { return &sigma2_;}

      //setter
      /** change the value of lambda_ */
      inline void setLambda(STK::Real const& lambda) {lambda_ = lambda;}
      /** change the value of sigma2_ */
      inline void setSigma2(STK::Real const& sigma2) {sigma2_ = sigma2;}

      //methods
      /**
       * update sigma2 and the lasso penalty
       * @param beta current estimates
       * @param normResidual ||y-X*beta||_2^2
       */
      void update(STK::CVectorX const& beta, STK::Real const& normResidual);

      /**
       * update the lasso penalty (for fixed sigma)
       * @param beta current estimates
       */
      void update(STK::CVectorX const& beta);

      /**
       * @param x a vector of length p_
       * @return the product invPenalty_*x
       */
      STK::CVectorX multInvPenalty(STK::CVectorX const& x) const;

      /**
       * @param x a vector of length p_
       * @return the product invPenalty_.sqrt()*x
       */
      STK::CVectorX multSqrtInvPenalty(STK::CVectorX const& x) const;

      /** penalty term
       *  @param beta current estimates
       *  @return t(beta) * penalty * beta
       */
      STK::Real penaltyTerm(STK::CVectorX const& beta) const;


    protected:
      /** update the penalty
       *  @param beta current estimates
       */
      void updatePenalty(STK::CVectorX const& beta);

      /** update sigma2
       *  @param beta current estimates
       *  @param normResidual ||y-X*beta||_2^2
       */
      void updateSigma2(STK::CVectorX const& beta, STK::Real const& normResidual);

    private:
      ///value associated to the l1 penalty of estimates
      STK::Real lambda_;
      ///diag(E[1/tau_i^2])^-1 =diag(|beta|/lambda_)
      STK::Array2DDiagonal<STK::Real> sqrtInvPenalty_;
      ///variance
      STK::Real sigma2_;
//      ///number of sample
//      int n_;
//      ///number of covariates
//      int p_;
  };
}

#endif /* LASSOPENALTY_H_ */
