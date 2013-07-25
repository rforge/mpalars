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
 * created on: 3 avr. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file Fusion.h
 *  @brief In this file, we define the @c Fusion class.
 **/


#ifndef FUSION_H_
#define FUSION_H_

#include "../lars/Lars.h"

namespace HD
{
/**
 * This class transform the input data in order to use the lars algorithm for the fusion problem
 * (l1 norm penalization in the difference of successive coefficient).
 */
  class Fusion
  {
    public:
      //constructors
      /**
       * constructor
       * @param X matrix of size n*p, a row contains the values of each covariate for an individual.
       * @param y vector of length n containing the response
       */
      Fusion(STK::CArrayXX const& X, STK::CVectorX const& y);

      /**
       *
       * @param X matrix of size n*p, a row contains the values of each covariate for an individual.
       * @param y vector of length n containing the response
       * @param maxSteps number of maximum step to do
       * @param eps epsilon (for 0)
       * @param verbose if TRUE print some details
       */
      Fusion(STK::CArrayXX const& X, STK::CVectorX const& y, int maxSteps, STK::Real eps =STK::Arithmetic<STK::Real>::epsilon(), bool verbose = false);


      //getters
      /**@return path of the coefficients*/
      inline Path path() const {return path_;}

      /**@return Number of step of the algorithm*/
      inline int step() const {return step_;}

      /**
       * @param i step
       * @return the Pathstate object : the state of the path at the step i
       */
      inline PathState const& path(int i) const {return path_.states(i);}

      /**
       * @param i index of the step
       * @param j index of the coefficients
       * @return the value of the j-th coefficient at the step i
       */
      inline STK::Real coefficient(int i,int j) const {return path_.varCoeff(i,j);}

      /**
       * @param i index of the step
       * @param j index of the coefficients
       * @return the value of the j-th coefficient at the step i
       */
      inline int varIdx(int i,int j) const {return path_.varIdx(i,j);}

      /**
       * @param i index of the step
       * @return the value of l1norm at the i-th step
       */
      inline STK::Real l1norm(int i) const {return path_.l1norm(i);}

      /**@return the value of lambda */
      inline std::vector<STK::Real> lambda() const {return path_.lambda();}

      /** @return the historic of add and drop variable*/
      inline std::vector< std::pair<int,int> > evolution() const {return path_.evolution();}

      /** @return the intercept of the solution*/
      inline STK::Real mu() const {return mu_;}

      /** @return the ignored variable*/
      inline std::vector<bool> toIgnore() const {return toIgnore_;}

      /**
       * run the lars algorithm for solving the fusion problem on Z=X*L^-1 (L^-1 = lower triangular matrix of 1)
       */
      void run();

    protected:
      /**
       * change X in Z=X*L^-1 (L^-1 = lower triangular matrix of 1)
       */
      void computeZ();

    private:
      ///matrix of size n*p, a col = a covariate
      STK::CArrayXX X_;
      ///vector size n, response
      STK::CVectorX y_;
      ///maximum number of steps for the lars algorithm
      int maxSteps_;
      ///numerical zero
      STK::Real eps_;
      ///if TRUE, print some details
      bool verbose_;
      ///mean
      STK::Real mu_;
      ///number of step done by the lars algorithm
      int step_;
      ///solution path of the lars
      Path path_;
      ///ignored variables (due to correlation)
      std::vector<bool> toIgnore_;

  };
}


#endif /* FUSION_H_ */
