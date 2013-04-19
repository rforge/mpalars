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
 *  @brief In this file .
 **/


#ifndef FUSION_H_
#define FUSION_H_

#include "./Lars.h"

namespace MPA
{
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
      Fusion(STK::CArrayXX const& X, STK::CVectorX const& y, int maxSteps, Real eps =STK::Arithmetic<Real>::epsilon(), bool verbose = false);


      //getters
      /**@return path of the coefficients*/
      inline Path path() const {return path_;}

      /**@return Number of step of the algorithm*/
      inline int step() const {return step_;}
      inline Real coefficient(int i,int j) const {return path_.varCoeff(i,j);}
      inline int varIdx(int i,int j) const {return path_.varIdx(i,j);}
      inline Real lambda(int i) const {return path_.lambda(i);}
      inline PathState const& path(int i) const {return path_.states(i);}
      inline std::vector< std::pair<int,int> > evolution() const {return path_.evolution();}
      inline Real mu() const {return mu_;}
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
      STK::CArrayXX X_;//matrix of size n*p, a col = a covariate
      STK::CVectorX y_;//vector size n, response

      int maxSteps_;
      Real eps_;
      bool verbose_;

      int step_;
      Real mu_;
      std::vector<bool> toIgnore_;
      Path path_;


  };
}


#endif /* FUSION_H_ */
