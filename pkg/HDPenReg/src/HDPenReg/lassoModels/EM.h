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
 * created on: 30 mai 2013
 * Author:   Quentin Grimonprez
 **/

/** @file EM.h
 *  @brief In this file, definition of the EM class .
 **/


#ifndef EM_H_
#define EM_H_

#include "IAlgo.h"


namespace HD
{

  /**
   * This class runs an EM algorithm on a @c PenalizedModels object.
   * The stopping criterion is the convergence of the completed loglikelihood or the number of iterations.
   */
  class EM : public IAlgo
  {
    public:
      /**
       * Constructor
       * @param maxStep maximal number of steps of the algorithm
       * @param eps threshold for convergence of the completed loglikelihood
       */
      EM(int maxStep, STK::Real eps=STK::Arithmetic<STK::Real>::epsilon()) : IAlgo(), maxStep_(maxStep), eps_(eps) {}

      /**@return the number of step of the algorithm*/
      inline int step() const {return step_;};

      /**
       * run the EM algorithm on a PenalizedModels object
       * @param model pointer to a PenalizedModels object
       */
      bool run(PenalizedModels*& model)
      {
        try
        {
          //initialization
          step_ = 0;
          STK::Real diff = eps_+1;//for make at least one loop
          STK::Real llc = model ->llc(), llcOld(0);

          //we stop the algorithm after convergence of the completed loglikelihood or after reaching the maw number of step
          while( diff>eps_ && step_<maxStep_)
          {
            //stk_cout<<" Step  "<<step_<<std::endl;
            model -> eStep();
            model -> mStep();

            //difference between the loglikelihood of 2 successive steps.
            llcOld=llc;
            llc = model -> llc();
            diff = std::abs(llc-llcOld);

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

      /** @return the last error message**/
      inline STK::String const& error() const { return msg_error_;}

    private:
      ///maximum number of steps of the EM
      int maxStep_;
      ///step of the algorithm
      int step_;
      ///threshold for convergence of the complete loglikelihood
      STK::Real eps_;
      ///last error message
      STK::String msg_error_;
  };

}

#endif /* EM_H_ */
