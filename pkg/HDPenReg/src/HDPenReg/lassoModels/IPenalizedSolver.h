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
 * created on: 14 juin 2013
 * Author:   Quentin Grimonprez
 **/

/** @file IPenalizedSolver.h
 *  @brief In this file definition of the interface class @c IPenalizedSolver.
 **/


#ifndef IPENALIZEDSOLVER_H_
#define IPENALIZEDSOLVER_H_

#include "../../stkpp/include/STKpp.h"

namespace HD
{
  /** @ingroup lassoModels
   *  @brief The class IPenalizedSolver is an interface for the solver of the @c PenalizedModels M-step
   */
  class IPenalizedSolver
  {
    public:
      /**
       * Constructor
       * @param p_currentData pointer to the current data
       * @param p_currentSet pointer to the current set of non zero estimates
       */
      IPenalizedSolver(STK::CArrayXX const* p_currentData = 0, STK::Array2DVector<int> const* p_currentSet = 0)
      : p_currentData_(p_currentData)
      , p_currentSet_(p_currentSet)
      {};

      /**destructor*/
      virtual ~IPenalizedSolver() {} ;

      /**
       * run the solver
       * @return the new estimated value of beta
       */
      virtual STK::CVectorX run() = 0;

      //setter
      /**set the pointer to the current data (data reduce to covariates from current set)*/
      inline void setCurrentData(STK::CArrayXX const* p_currentData){p_currentData_=p_currentData;};
      /**set the pointer to the current set*/
      inline void setCurrentSet(STK::CArrayXX const* p_currentSet){p_currentData_=p_currentSet;};

    protected:
      ///const pointer on the current data (data reduce to variable of non-zero estimates)
      STK::CArrayXX const* p_currentData_;
      ///const pointer on the current set of variable (contains the index of non-zero estimates)
      STK::Array2DVector<int> const* p_currentSet_;

  };
}

#endif /* IPENALIZEDSOLVER_H_ */
