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
 * created on: 4 nov. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file LogisticLassoSolver.h
 *  @brief In this file .
 **/


#ifndef LOGISTICLASSOSOLVER_H_
#define LOGISTICLASSOSOLVER_H_

#include "IPenalizedSolver.h"
#include "LassoPenalty.h"
#include "LogisticLassoPenalty.h"

namespace HD
{

  /** @ingroup lassoModels
   *  @brief This class inherits from the @c IPenalizedSolver class. It implements the way to solve the Mstep for a logistic lasso penalty
   */
  class LogisticLassoSolver : public IPenalizedSolver
  {
    public:

      /**default constructor*/
      LogisticLassoSolver();

      /**
       * Constructor
       * @param p_data pointer to the current Data
       * @param beta initial solution
       * @param p_y pointer to the response
       * @param threshold threshold for shrinkage
       * @param p_solver pointer to the solver
       * @param p_penalty pointer to the lasso penalty
       */
      LogisticLassoSolver(STK::CArrayXX const* p_data, STK::CVectorX const& beta, STK::CVectorX const* p_y = 0,  STK::Real const& threshold = 1e-10,
                  STK::CG<LassoMultiplicator,STK::CVectorX,InitFunctor>* p_solver = 0, LogisticLassoPenalty* p_penalty = 0 );

      /**destructor*/
      virtual ~LogisticLassoSolver() {};

      /**Solve the M-step with a conjugate gradient
       * @return the completed loglikelihood
       * */
      STK::Real run(bool const& burn = true);

      /**run the update of the penalty*/
      void update();

      /**Initialization of the solver*/
      void initializeSolver();

      //getter
      /**@return the pointer to the penalty*/
      inline LogisticLassoPenalty* p_penalty() const { return p_penalty_;}
      /**@return a pointer to the CG solver*/
      inline STK::CG<LassoMultiplicator,STK::CVectorX,InitFunctor>* p_solver() {return p_solver_;}

      //setter
      /** set the conjugate gradient solver
       * @param p_solver pointer to the solver
       */
      inline void setSolver(STK::CG<LassoMultiplicator,STK::CVectorX,InitFunctor>* p_solver) {p_solver_ = p_solver;}

      /**
       * set the LassoPenalty
       * @param p_penalty pointer to the penakty
       */
      inline void setPenalty(LogisticLassoPenalty* p_penalty) {p_penalty_ = p_penalty;}
      /**
       * set the threshold
       * @param threshold threshold for shrinkage to 0
       */
      inline void setThreshold(STK::Real threshold) {threshold_ = threshold;}


    protected:
      /**Thresholding of the new estimates : estimated coefficients < threshold_ become 0*/
      void thresholding();

      /** update all the current variables*/
      void updateCurrent();

      /**Update the currentBeta_ and currentData_*/
      void updateCurrentData();

      /** Computation of the completed loglikelihood*/
      STK::Real computeLlc();

    private:
      ///pointer to the conjugate gradient with LassoMultiplicator
      STK::CG<LassoMultiplicator,STK::CVectorX,InitFunctor>* p_solver_;
      ///t(X) * y
      STK::CVectorX Xty_;
      ///b from ax=b for CG
      STK::CVectorX b_;
      ///number of active variables in the current set
      int nbActiveVariables_;
      ///threshold under we consider a beta equal to 0
      STK::Real threshold_;
      ///pointer to the lasso penalty
      LogisticLassoPenalty* p_penalty_;

  };
}


#endif /* LOGISTICLASSOSOLVER_H_ */
