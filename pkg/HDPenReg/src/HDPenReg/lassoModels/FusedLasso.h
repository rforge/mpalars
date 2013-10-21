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
 * created on: 2 oct. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file FusedLasso.h
 *  @brief In this file, definition of the class @c FusedLasso .
 **/


#ifndef FUSEDLASSO_H_
#define FUSEDLASSO_H_

#include "PenalizedModels.h"
#include "FusedLassoPenalty.h"
#include "FusedLassoSolver.h"

namespace HD
{
  class FusedLasso;

  template<>
  struct ModelTraits<FusedLasso>
  {
    typedef FusedLassoSolver Solver;
    typedef FusedLassoPenalty Penalty;
    typedef FusedLassoMultiplicator Multiplicator;
    typedef STK::CG<FusedLassoMultiplicator,STK::CVectorX, InitFunctor> CG;
  };

  /**
   * Class FusedLasso derived from @c PenalizedModels.
   * This class constructs a FusedLasso Model to be solved by an @c EM algorithm.
   */
  class FusedLasso : public PenalizedModels<FusedLasso>
  {
    public:
    /**default constructor*/
      FusedLasso()
      : PenalizedModels<FusedLasso>()
      {
        STK::CVectorX Xty(0);
        // creation of the fusedlasso penalty
        FusedLassoPenalty* p_penalty = new FusedLassoPenalty();
        p_penalty_ = p_penalty;
        //creation functor for CG
        FusedLassoMultiplicator mult(0, p_penalty_->p_mainDiagonal(), p_penalty_->p_offDiagonal(), p_penalty_->p_sigma2());

        mult_ = mult;
        //create CG
        STK::CG<FusedLassoMultiplicator,STK::CVectorX, InitFunctor>* p_gcsolver = new STK::CG<FusedLassoMultiplicator,STK::CVectorX,InitFunctor>(mult_,Xty);

        p_gcsolver_ = p_gcsolver;

        //create solver for fused lasso
        FusedLassoSolver* p_fusedlassosolver = new FusedLassoSolver;
        p_fusedlassosolver->setPenalty(p_penalty_);
        p_fusedlassosolver->setSolver(p_gcsolver_);

        //add solver to lasso
        p_solver_ = p_fusedlassosolver;

        mult_.p_data_ = p_currentData();
      };

      /**
       * constructor
       * @param p_data pointer to the data
       * @param p_y pointer to the response
       * @param lambda1 value of parameter associated to the l1 penalty
       * @param lambda2 value of parameter associated to the l1 penalty of successive coefficients
       * @param threshold threshold for setting coefficient to 0
       * @param epsCG epsilon for CG convergence
       */
      FusedLasso(STK::CArrayXX const* p_data, STK::CVectorX const* p_y, STK::Real lambda1, STK::Real lambda2, STK::Real threshold, STK::Real epsCG)
      : PenalizedModels<FusedLasso>(p_data, p_y)
      {
        STK::CVectorX Xty(nbVar());
        Xty = p_data_->transpose() * *p_y;

        STK::CVectorX beta0(Xty.sizeRows());
        for(int i = 1; i <= Xty.sizeRows();i++)
          beta0[i] = Xty[i]/((*p_data).col(i).norm2());

        beta_ = beta0;

        // creation fusedlasso penalty
        FusedLassoPenalty* p_penalty = new FusedLassoPenalty(lambda1, lambda2, threshold);
        p_penalty_ = p_penalty;
        //creation functor for CG
        FusedLassoMultiplicator mult(0, p_penalty_->p_mainDiagonal(), p_penalty_->p_offDiagonal(), p_penalty_->p_sigma2());

        mult_ = mult;
        //create CG
        STK::CG<FusedLassoMultiplicator,STK::CVectorX, InitFunctor>* p_gcsolver = new STK::CG<FusedLassoMultiplicator,STK::CVectorX,InitFunctor>(mult_, Xty,0 ,epsCG);

        p_gcsolver_=p_gcsolver;
        //create solver for fused lasso
        FusedLassoSolver* p_fusedlassosolver = new FusedLassoSolver(p_data_, beta_, p_y_, p_gcsolver_, p_penalty_, threshold);

        //add the penalty
        p_fusedlassosolver->setPenalty(p_penalty_);

        //add solver to lasso
        p_solver_ = p_fusedlassosolver;

        InitFunctor init(p_solver_->p_currentBeta());
        init_ = init;
        (p_solver_->p_solver())->setInitFunctor(&init_);

        mult_.p_data_ = p_currentData();
      }

      /** destructor*/
      ~FusedLasso()
      {
        if(p_gcsolver_) delete p_gcsolver_;
      }

      /**
       * set the lasso regularization parameter
       * @param lambda1
       */
      inline void setLambda1(STK::Real const& lambda1) {(p_solver_->p_penalty())->setLambda1(lambda1);}
      /**
       * set the fusion regularization parameter
       * @param lambda2
       */
      inline void setLambda2(STK::Real const& lambda2) {(p_solver_->p_penalty())->setLambda2(lambda2);}

      /**
       * set the threshold for segments
       * @param threshold
       */
      inline void setThreshold(STK::Real const& threshold) {p_solver_->setEps(threshold);}

      /**
       * set the epsilon for avoid 0 to denominator
       * @param epsilon
       */
      inline void setEps(STK::Real const& eps) {(p_solver_->p_penalty())->setEps(eps);}

      /**
       * set the epsilon for the CG
       * @param epsilon epsilon for the convergence of CG
       */
      inline void setCGEps(STK::Real const& eps) {p_gcsolver_->setEps(eps);}

      /**initialize the containers of all subclasses*/
      void initializeModel()
      {
        //computation of the initial solution : beta0
        STK::CVectorX Xty(nbVar());
        Xty = p_data_->transpose() * (*p_y_);

        STK::CVectorX beta0(Xty.sizeRows());
        for(int i = 1; i <= Xty.sizeRows();i++)
          beta0[i] = Xty[i]/((*p_data_).col(i).norm2());

        beta_ = beta0;

        //set the parameters of the FusedLassoSolver
        p_solver_->setData(p_data_);
        p_solver_->setY(p_y_);
        p_solver_->setBeta(beta_);
        p_solver_->initializeSolver();

        //initialization functor for the CG
        InitFunctor init(p_solver_->p_currentBeta());
        init_ = init;
        p_gcsolver_->setInitFunctor(&init_);

      }

    private:
      /// multiplicator for conjugate gradient
      FusedLassoMultiplicator mult_;
      /// conjugate gradient for solver
      STK::CG<FusedLassoMultiplicator,STK::CVectorX, InitFunctor>* p_gcsolver_;
      /// initial functor for CG
      InitFunctor init_;
  };
}

#endif /* FUSEDLASSO_H_ */
