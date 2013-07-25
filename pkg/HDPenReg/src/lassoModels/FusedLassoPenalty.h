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

    Contact : Serge.Iovleff@stkpp.org
*/

/*
 * Project:  MPAGenomics::
 * created on: 21 juin 2013
 * Author:   Quentin Grimonprez
 **/

/** @file FusedLassoPenalty.h
 *  @brief In this file .
 **/


#ifndef FUSEDLASSOPENALTY_H_
#define FUSEDLASSOPENALTY_H_

#include "PenalizedModels.h"
#include "IPenalty.h"

namespace HD
{
  struct FusedLassoMultiplicator
  {
      STK::CVectorX operator()(STK::CVectorX &x) const
      {
        STK::Array2DDiagonal<STK::Real> sigIdx(p_data_->sizeColsImpl(),*p_sigma2_);
        STK::CVectorX a(x.sizeRowsImpl());

        for(int i = 1; i <= x.sizeRowsImpl(); i++)
        {
          sigIdx[i] *= x[i];
          a[i] = std::sqrt((*p_invPenalty_)[i]) * x[i];//invD*x
        }

        //a = (p_invPenalty->sqrt()) * x;
        a = (*p_data_) * a;//X*invD*x
        a = p_data_->transpose() * a;//tX*X*invD*x
        //a = (p_invPenalty->sqrt()) * a;
        for(int i = 1; i <= x.sizeRowsImpl(); i++)
          a[i] = std::sqrt((*p_invPenalty_)[i]) * a[i];//invD*tX*X*invD*x

        a += sigIdx;//sigI*x+invD*tX*X*invD*x

        return   a ;
      }
      STK::CArrayXX A()
      {
        STK::Array2DDiagonal<STK::Real> sigId(p_data_->sizeColsImpl(),*p_sigma2_),sqrtInvPen(p_data_->sizeColsImpl());

        for(int i = 1; i< p_data_->sizeColsImpl();i++)
          sqrtInvPen[i] = std::sqrt((*p_invPenalty_)[i]);
        STK::CArrayXX b;
        b =  *p_invPenalty_ * ( p_data_->transpose() * (*p_data_) ) + sigId;

        return b;
      }
      FusedLassoMultiplicator(STK::CArrayXX const* p_data,STK::Array2DDiagonal<STK::Real> const* p_invPenalty,STK::Real const* p_sigma2)
                        : p_data_(p_data), p_invPenalty_(p_invPenalty), p_sigma2_(p_sigma2)
      {
      }

      STK::CArrayXX const* p_data_;
      STK::Array2DDiagonal<STK::Real> const* p_invPenalty_;
      STK::Real const* p_sigma2_;
  };



  class FusedLassoPenalty : public IPenalty
  {
    public:
      /** Constructor
       *  @param lambda penalization parameter for the l1-norm of the estimates
       *  @param n size of sample
       *  @param p size of Penalty (number of covariates)
       */
      FusedLassoPenalty(STK::Real lambda, int n, int p);

      /** Copy constructor
       *  @param penalty LassoPenalty object to copy
       */
      FusedLassoPenalty(LassoPenalty const& penalty);

      /** destructor */
      inline virtual ~LassoPenalty() {};

      /**clone*/
      FusedLassoPenalty* clone() const;

      //getter
      /**@return lambda parameter of the lasso */
      inline STK::Real const& lamba() const {return lambda_;}
      /**@return invPenalty diagonal matrix containing |beta_i| / lambda */
      inline STK::Array2DDiagonal<STK::Real> const& invPenalty() const {return invPenalty_;}
      /**@return n size of sample */
      inline int const& n() const {return n_;}
      /**@return p number of covariates */
      inline int const& p() const {return p_;}
      /**@return sigma2 variance of the response*/
      inline STK::Real const& sigma2() const { return sigma2_;}

      inline STK::Array2DDiagonal<STK::Real> const* p_invPenalty() const {return &invPenalty_;}
      inline STK::Real const*  p_sigma2() const { return &sigma2_;}

      STK::CArrayXX A() const;

      /**
       * update sigma2 and the lasso penalty
       * @param beta current estimates
       * @param normResidual ||y-X*beta||_2^2
       */
      void update(STK::CVectorX const& beta, STK::Real const& normResidual);

      STK::CVectorX multInvPenalty(STK::CVectorX const& x) const;
      STK::CVectorX multSqrtInvPenalty(STK::CVectorX const& x) const;


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
      STK::Real lambda_;
      STK::Array2DDiagonal<STK::Real> invPenalty_;//diag(E[1/tau_i^2])^-1
      STK::Real sigma2_;
      int n_;
      int p_;
  };
}

#endif /* FUSEDLASSOPENALTY_H_ */
