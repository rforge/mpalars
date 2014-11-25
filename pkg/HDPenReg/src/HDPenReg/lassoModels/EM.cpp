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
 * created on: 17 december 2013
 * Author:   Quentin Grimonprez
 **/

/** @file EM.cpp
 *  @brief In this file, methods of class @c EM are implemented.
 **/



#include "EM.h"

namespace HD
{
  /*
   * Constructor
   * @param maxStep maximal number of steps of the algorithm
   * @param eps threshold for convergence of the completed loglikelihood
   */
  EM::EM(int maxStep, int burn, STK::Real eps) :
      IAlgo(),
      burn_(burn)
  {
    eps_ = eps;
    maxStep_ = maxStep;
  }


  /**
   * run the EM algorithm on a PenalizedModels object
   * @param model pointer to a PenalizedModels object
   */
  bool EM::run(PenalizedModelsBase* model)
  {
    try
    {
      //initialization
      step_ = 0;
      STK::Real diff = eps_+1;//for make at least one loop
      STK::Real llc = std::numeric_limits<STK::Real>::max(), llcOld(0);

      //we stop the algorithm after convergence of the completed loglikelihood or after reaching the max number of step
      while( diff>eps_ && step_<maxStep_)
      {
        model -> eStep();
        model -> mStep( (burn_<=step_) );

        //difference between the loglikelihood of 2 successive steps.
        llcOld=llc;
        llc = model -> lnLikelihood();
        diff = std::abs((llc-llcOld)/llcOld);

        step_++;
      }

    }
    catch(const STK::Exception& e)
    {
      msg_error_ = e.error();
      return false;
    }

    return true;
  }

}
