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
 * created on: 14 juin 2013
 * Author:   Quentin Grimonprez
 **/

/** @file IPenalizedSolver.h
 *  @brief In this file definition of the interface class @c IPenalizedSolver.
 **/


#ifndef IPENALIZEDSOLVER_H_
#define IPENALIZEDSOLVER_H_

#include "../../stkpp/include/STKpp.h"
#include "IPenalty.h"


namespace HD
{

  /**
   * Functor for initialization in the conjugate gradient
   */
  struct InitFunctor
  {
      /**
       * Operator
       * @return pointer on the value for initialization
       */
      STK::CVectorX operator()() const
      { return *p_x_;}

      /**
       * Constructor
       * @param p_x pointer on the value for initialization
       */
      InitFunctor(STK::CVectorX const* p_x = 0) : p_x_(p_x) {};

      ///pointer on the value for initialization
      STK::CVectorX const* p_x_;
  };

  /** @ingroup lassoModels
   *  @brief The class IPenalizedSolver is an interface for the solver of the @c PenalizedModels M-step
   */
  class IPenalizedSolver
  {
    public:
      /**default constructor*/
      IPenalizedSolver()
      : currentData_()
      , currentBeta_()
      , currentSet_()
      , p_data_(0)
      , p_y_(0)
      , beta_()
      {
      }

      /**
       * Constructor
       * @param beta value for initializing beta
       * @param p_data pointer to the data
       * @param p_y pointer to the response
       */
      IPenalizedSolver(STK::CVectorX const& beta, STK::CArrayXX const* p_data = 0, STK::CVectorX const* p_y = 0)
                           : currentData_(*p_data)
                           , currentBeta_(beta)
                           , currentSet_(beta.sizeRows())
                           , p_data_(p_data)
                           , p_y_(p_y)
                           , beta_(beta)
       {
         for(int i = 1; i <= currentSet_.sizeRows(); i++)
           currentSet_[i] = i;
       };

       /**destructor*/
       virtual ~IPenalizedSolver() {} ;

       /**
        * run the solver (Mstep)
        * @return the new estimated value of beta
        */
       virtual STK::Real run(bool const& burn = true) = 0;
       /**run the update of the penalty (Estep)*/
       virtual void update() = 0;
       /** initialize all the containers of the class */
       virtual void initializeSolver() = 0;

       //setter
       /**set the pointer to the current data (data reduce to covariates from current set)*/
       inline void setData(STK::CArrayXX const* p_data) {p_data_ = p_data;};
       /**set the pointer to the response*/
       inline void setY(STK::CVectorX const* p_y) {p_y_ = p_y;};
       /**set the pointer to the current beta*/
       inline void setBeta(STK::CVectorX const& beta)
       {
         beta_ = beta;
         currentBeta_ = beta_;
       };


       //getter
       /**@return beta_*/
       inline STK::CVectorX  const& beta() const {return beta_;};
       /**@return currentData_*/
       inline STK::CArrayXX const& currentData() const {return currentData_;};
       /**@return currentBeta_*/
       inline STK::CVectorX const& currentBeta() const {return currentBeta_;};
       /**@return currentData_*/
       inline STK::CArrayXX const* p_currentData() const {return &currentData_;};
       /**@return currentBeta_*/
       inline STK::CVectorX const* p_currentBeta() const {return &currentBeta_;};
       ///@return currentSet_
       inline STK::Array2DVector<int> const& currentSet() const {return currentSet_;};


    protected:
       /** compute the loglikelihood
        * @return the loglikelihood of the current step
        */
       virtual STK::Real computeLlc() = 0;

    protected:
       ///current data
       STK::CArrayXX currentData_;
       ///Current beta
       STK::CVectorX currentBeta_;
       ///current set of variable
       STK::Array2DVector<int> currentSet_;
       ///const pointer on the data
       STK::CArrayXX const* p_data_;
       ///const pointer on the response
       STK::CVectorX const* p_y_;
       ///const pointer on beta_
       STK::CVectorX beta_;

   };


}

#endif /* IPENALIZEDSOLVER_H_ */
