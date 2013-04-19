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
 * created on: 5 févr. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file PathState.h
 *  @brief In this file, we define the final class @c PathState.
 **/


#ifndef PATHSTATE_H_
#define PATHSTATE_H_

#include "stkpp/include/STKpp.h"
#include <utility>

namespace MPA
{
  class PathState
  {
    public:
      //constructors
      /** default constructor*/
      PathState();

      /**
       * constructor with reserve size
       * @param nbMaxVariable maximal number of variable to potentially stock
       */
      PathState(int nbMaxVariable);

      //setters
      /**set indexVariables_ to a new value
       * @param indexVariables the new value of indexVariables_
       * */
      inline void setCoefficients(STK::Array2DVector< std::pair<int,Real> > const& coefficients){coefficients_=coefficients;}
      inline void setLambda(Real const& lambda) {lambda_=lambda;}

      //getters
      /**@return lambda_*/
      inline Real const lambda() const {return lambda_;}
      /**@return coefficients_*/
      inline STK::Array2DVector< std::pair<int,Real> > const& coefficients() const {return coefficients_;}
      /**@return coefficients_[i]*/
      inline std::pair<int,Real> const& coefficients(int i) const {return coefficients_[i];}
      /**@return coefficients_[i].first*/
      inline int  varIdx(int i) const {return coefficients_[i].first;}
      /**@return coefficients_[i].first*/
      inline Real varCoeff(int i) const {return coefficients_[i].second;}
      /**@return size of the vector*/
      inline int size() const {return coefficients_.size();}

      //methods
      /** update coefficients with values specified in parameters
       * @param indexVariables vector containing the index of active variables
       * @param coefficients vector containing the value of estimates for active variables
       */
      void update(STK::Array2DVector<int> const& indexVariables,STK::Array2DVector<Real> const& coefficients);

      /**
       * update of the coefficients of the previous step
       * @param w direction of the update
       * @param gamma step of the update
       */
      void update(STK::Array2DVector<Real> const& w, Real gamma);

      /**
       * update of the coefficients of the previous state with a variable to drop and a variable to add
       * @param w direction of the update
       * @param gamma step of the update
       * @param addIdxVar index of the variable to add
       * @param dropIdx index (in the vector of coefficients of the previous step) of the variable to delete
       */
      void addWithDropUpdate(STK::Array2DVector<Real> const& w, Real gamma, int addIdxVar, int dropIdx);

      /**
       * update of the coefficients of the previous state with a variable to drop
       * @param w direction of the update
       * @param gamma step of the update
       * @param dropIdx index (in the vector of coefficients of the previous step) of the variable to delete
       */
      void dropAfterDropUpdate(STK::Array2DVector<Real> const& w, Real gamma, int dropIdx);

      /**
       * update of the coefficients of the previous state with a new variable
       * @param w direction of the update
       * @param gamma step of the update
       * @param addIdxVar index of the variable to add
       */
      void addUpdate(STK::Array2DVector<Real> const& w, Real gamma, int addIdxVar);

      /**print coefficients*/
      void printCoeff() const;


    private:
      //attributes
      STK::Array2DVector< std::pair<int,Real> > coefficients_;
      Real lambda_;
  };
}

#endif /* PATHSTATE_H_ */