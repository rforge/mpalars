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

/** @file EnetSolver.h
 *  @brief In this file, definition of the @c EnetSolver class.
 **/


#ifndef ENETSOLVER_H_
#define ENETSOLVER_H_

#include "IPenalizedSolver.h"
#include "EnetPenalty.h"

namespace HD
{
  /**
   * This class inherits from the @c IPenalizedSolver class.
   * It implements the way to solve the Mstep of an elastic net problem.
   */
  class EnetSolver : public IPenalizedSolver
  {
    public:
      /**
       * Constructor
       * @param p_currentData pointer to the current Data
       * @param p_currentSet pointer to the current Set
       * @param Xty t(X) * y
       * @param p_solver pointer to the solver
       * @param p_penalty pointer to the lasso penalty
       */
      EnetSolver(STK::CArrayXX const* p_currentData, STK::Array2DVector<int> const* p_currentSet, STK::CVectorX const& Xty, STK::CG<EnetMultiplicator,STK::CVectorX,InitFunctor>* p_solver = 0, EnetPenalty* p_penalty = 0 )
                : IPenalizedSolver()
                , p_solver_(p_solver)
                , p_penalty_(p_penalty)
                , Xty_(Xty)
      {
      }
      virtual ~EnetSolver() {};

      /**
       * Solve the M-step with a conjugate gradient
       * @return new estimate of beta
       */
      STK::CVectorX run()
      {
        //compute the b of the linear system Ax=b
        STK::CVectorX XtyTemp(p_currentSet_->sizeRows());
        for(int i = 1; i <= XtyTemp.sizeRows(); i++)
          XtyTemp[i]=Xty_[(*p_currentSet_)[i]];

        //set the b of the linear system Ax=b
        p_solver_->setB(p_penalty_->multSqrtInvPenalty(Xty_));

        //run the conjugate gradient
        bool convergence;
        convergence=p_solver_->run();

        //backtransform the solution x to beta
        return p_penalty_->multSqrtInvPenalty(p_solver_->x());
      }


      /** set the conjugate gradient solver*/
      inline void setSolver(STK::CG<EnetMultiplicator,STK::CVectorX,InitFunctor>* p_solver) {p_solver_=p_solver;}
      /** set the pointer to the lassopenalty*/
      inline void setPenalty(EnetPenalty* p_penalty) {p_penalty_=p_penalty;}
      /**set Xty_*/
      inline void setXty(STK::CVectorX const& Xty) {Xty_=Xty;}


    private:
      ///pointer to the conjugate gradient with LassoMultiplicator
      STK::CG<EnetMultiplicator,STK::CVectorX,InitFunctor>* p_solver_;
      ///pointer to the Enet penalty
      EnetPenalty* p_penalty_;
      ///t(X) * y
      STK::CVectorX Xty_;
  };
}

#endif /* ENETSOLVER_H_ */
