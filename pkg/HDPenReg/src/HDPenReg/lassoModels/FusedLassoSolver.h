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
 * created on: 17 sept. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file FusedLassoSolver.h
 *  @brief In this file, definition of the @c FusedLassoSolver class .
 **/


#ifndef FUSEDLASSOSOLVER_H_
#define FUSEDLASSOSOLVER_H_

#include "IPenalizedSolver.h"
#include "FusedLassoPenalty.h"

namespace HD
{


  /** @ingroup lassoModels
   *  @brief This class inherits from the @c IPenalizedSolver class. It implements the way to solve the Mstep for a Fused lasso penalty
   */
  class FusedLassoSolver : public IPenalizedSolver
  {
    public:

      /**default constructor*/
      FusedLassoSolver();
      /**
       * Constructor
       * @param p_data pointer to the full data
       * @param beta initial solution of the problem
       * @param p_y pointer to the response associated with the data
       * @param p_solver pointer to the solver
       * @param p_penalty pointer to the Fused lasso penalty
       * @param eps tolerance for zero
       */

      FusedLassoSolver(STK::CArrayXX const* p_data, STK::CVectorX const& beta, STK::CVectorX const* p_y = 0,
                       STK::CG<FusedLassoMultiplicator,STK::CVectorX,InitFunctor>* p_solver = 0, FusedLassoPenalty* p_penalty = 0,
                       STK::Real eps = 1e-8);


      /**destructor*/
      ~FusedLassoSolver() {};

      /**
       * Solve the M-step with a conjugate gradient
       * @return the completed loglikelihood
       */
      STK::Real run(bool const& burn = true);

      /**run the update of the penalty*/
      void update();

      /**initialization of the solver*/
      void initializeSolver();

      //setter and getter
      /** set the conjugate gradient solver
       * @param p_solver pointer to a conjugate gradient solver
       * */
      inline void setSolver(STK::CG<FusedLassoMultiplicator,STK::CVectorX,InitFunctor>* p_solver) {p_solver_ = p_solver;}
      /** set the pointer to the FusedLassopenalty
       * @param p_penalty pointer to the FusedLassopenalty
       * */
      inline void setPenalty(FusedLassoPenalty* p_penalty) {p_penalty_ = p_penalty;}
      /** set the tolerance for the difference of estimates
       * @param eps new tolerance
       */
      inline void setEps(STK::Real eps) {eps_ = eps;}

      /** @return segment_ */
      inline std::vector<STK::Range> segment() const {return segment_;}

      /**@return pointer to the penalty*/
      inline FusedLassoPenalty* p_penalty() {return p_penalty_;}
      /**@return pointer to the solver*/
      inline STK::CG<FusedLassoMultiplicator,STK::CVectorX,InitFunctor>* p_solver() {return p_solver_;}


    protected:
      /**
       * update currentSet_, segment_ and currentBeta_
       * @return a boolean, if true there is a changement in the currentSet
       */
      bool updateCurrent();

      /**update currentData_*/
      void updateCurrentData();

      /** update beta_*/
      void updateBeta();

      /**Computation of the completed log-likelihood
       * @return the current completed loglikelihood
       */
      STK::Real computeLlc();

    private:
      ///pointer to the conjugate gradient with LassoMultiplicator
      STK::CG<FusedLassoMultiplicator,STK::CVectorX,InitFunctor>* p_solver_;
      ///currentData_.transpose() * (*p_y_)
      STK::CVectorX currentXty_;
      ///vector containing the segment corresponding to every index of beta
      std::vector<STK::Range> segment_;
      /// number of sample in the data
      int n_;
      ///size of currentBeta_
      int nbActiveVariables_;
      ///eps_ tolerance for the difference of estimates
      STK::Real eps_;
      ///pointer to the penalty
      FusedLassoPenalty* p_penalty_;
  };
}


#endif /* FUSEDLASSOSOLVER_H_ */
