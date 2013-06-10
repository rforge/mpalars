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
 * created on: 4 juin 2013
 * Author:   Quentin Grimonprez
 **/

/** @file Cvlars.h
 *  @brief In this file .
 **/


#ifndef CVLARS_H_
#define CVLARS_H_

#include "Lars.h"
#include "functions.h"

namespace MPA
{
  class Cvlars
  {
    public:
      /**
       * Constructor with no index ( it will be a sequence from 0 to 1 by 0.01)
       * @param X matrix of data, a row=a individual
       * @param y response
       * @param k number of folds
       * @param maxSteps number of maximum step to do
       * @param eps epsilon (for 0)
       */
      Cvlars(STK::CArrayXX const& X, STK::CVectorX const& y, int k, int maxSteps, STK::Real eps = STK::Arithmetic<STK::Real>::epsilon());

      /**
       * Constructor
       * @param X matrix of data, a row=a individual
       * @param y response
       * @param k number of folds
       * @param index vector with real between 0 and 1 (ratio (norm coefficient)/max(norm coefficient) for which we compute the prediction error)
       * @param maxSteps number of maximum step to do
       * @param eps epsilon (for 0)
       */
      Cvlars(STK::CArrayXX const& X, STK::CVectorX const& y, int k, std::vector<double> const& index, int maxSteps, STK::Real eps = STK::Arithmetic<STK::Real>::epsilon());

      /**
       * run the cross validation
       */
      void run();

      //getter
      /** @return return the prediction error for each index*/
      inline STK::CVectorX const& cv() const {return cv_;}
      /** @return return the standard deviation of prediction error for each index*/
      inline STK::CVectorX const& cvError() const {return cvError_;}
      /** @return return the index*/
      inline std::vector<double> const& index() const {return index_;}

    private:
      /**
       * create a random partition in k folds
       */
      void partition();
      /**
       * run cross validation for folds from idxStartFold to idxEndFold
       * @param idxStartFold index of the first fold
       * @param idxEndFold index of the last fold
       */
      void subrun(int idxStartFold,int idxEndFold);

    private:
      STK::CArrayXX const* p_X_;
      STK::CVectorX const* p_y_;
      std::vector<int> partition_;
      std::vector<int> sizePartition_;
      std::vector<double> index_;
      STK::CArrayXX residuals_;
      STK::CVectorX cv_;
      STK::CVectorX cvError_;
      int k_;
      int n_;
      int p_;
      int maxSteps_;
      STK::Real eps_;
  };
}

#endif /* CVLARS_H_ */
