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
 * created on: 9 oct. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file CVLasso.h
 *  @brief In this file, definition of the class @c CVLasso derived from class @c CV .
 **/


#ifndef CVLASSO_H_
#define CVLASSO_H_

#include "CV.h"
#include "EM.h"
#include "Lasso.h"
#include "IMeasure.h"

namespace HD
{
/**
 * Class derived from @c CV, implementing the cross validation for @c Lasso with @c EM algorithm.
 * This class contains setters and implementation of the pure virtual runModel method from @c CV.
 */
  template<class LassoModel>
  class CVLasso : public CV
  {
    public:
      /**default constructor*/
      CVLasso() : CV(){};

      /**initialize containers and class*/
      void initialize() {initializeCV();};

      /**set the epsilon for the convergence of the @c EM algorithm*/
      inline void setEps(STK::Real const& eps) {eps_ = eps;}
      /**set the maximum number of step of the @c EM algorithm*/
      inline void setMaxStep(int const& maxStep) {maxStep_ = maxStep;}
      /**set the number of burn steps of the @c EM algorithm*/
      inline void setBurn(int const& burn) {maxStep_ = burn;}
      /**set the epsilon for the convergene of the conjugate gradient (@c CG)*/
      inline void setEpsCG(STK::Real const& epsCG) {epsCG_ = epsCG;}
      /**set the threshold of the @c LassoSolver*/
      inline void setThreshold(STK::Real const& threshold) {threshold_ = threshold;}
      /**set the type of measure for evaluate the model*/
      inline void setTypeMeasure(IMeasure* p_typeMeasure) {p_typeMeasure_ = p_typeMeasure;}


    protected:
      /**
       * run lasso on all values of index for the given control data
       * @param i index of the actual delete folder
       * @param XTest data test
       * @param yTest response test
       * @param p_XControl pointer to the data control
       * @param p_yControl pointer to the response control
       */
      void runModel(int i, STK::CArrayXX const& XTest, STK::CVectorX const& yTest, STK::CArrayXX const* p_XControl, STK::CVectorX const* p_yControl)
      {
        STK::CVectorX yPred(sizePartition_[i] );

        //create em algorithm
        EM algo(maxStep_,burn_,eps_);

        //create model
        LassoModel lasso;
        //set parameters of the model
        lasso.setCGEps(epsCG_);
        lasso.setThreshold(threshold_);
        lasso.setP_data(p_XControl);
        lasso.setP_y(p_yControl);

        //run the lasso on all value of index
        for(int s = 1 ; s <= (int) index_.size(); s++)
        {
          //set the new value of lambda to test
          lasso.setLambda(index_[s-1]);
          //initialize the model
          lasso.initializeModel();
          //run the algo
          algo.run(&lasso);

          //we compute the prediction of the y associated to XTest
          yPred = XTest * lasso.beta();
          //compute the residuals
          measure_(s,i+1) = p_typeMeasure_->measure(yTest,yPred);
        }
      }

    private:
      /// eps for EM algorithm convergence
      STK::Real eps_;
      /// threshold for set value to 0 in LassoPenalty
      STK::Real threshold_;
      /// eps for conjugate gradient convergence
      STK::Real epsCG_;
      /// maximum number of step for EM algorithm
      int maxStep_;
      /// burn for EM algorithm
      int burn_;
      ///type of measure
      IMeasure* p_typeMeasure_;
  };


}


#endif /* CVLASSO_H_ */
