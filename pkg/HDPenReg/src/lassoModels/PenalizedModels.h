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
 * created on: 23 mai 2013
 * Author:   Quentin Grimonprez
 **/

/** @file PenalizedModels.h
 *  @brief In this file, we find the definition of the PenalizedModels class.
 **/


#ifndef PENALIZEDMODELS_H_
#define PENALIZEDMODELS_H_

#include <stkpp/projects/StatModels/include/STK_ILatentModel.h>
#include <stkpp/include/STKpp.h>
#include "IPenalizedSolver.h"
#include "IPenalty.h"

namespace HD
{
  /** @ingroup lassoModels
   *  @brief The class PenalizedModels inherits from the template class @c ILatentModel.
   *  It contains implementations of Estep and Mstep for run an EM algorithm.
   */
  class PenalizedModels : public STK::ILatentModel<STK::CArrayXX>
  {
    public:
      /**
       * Constructor
       * @param data each row contains values of covariates for an individuals
       * @param y response of each individuals
       * @param beta Initialization of beta_ (least-Square solution)
       * @param p_penalty pointer on the penalty of the model
       * @param p_solver pointer on the solver to use in the Mstep
       * @param threshold tolerance for thresholding the solution
       *
       */
      PenalizedModels(STK::CArrayXX const* data, STK::CVectorX const& y, STK::CVectorX const& beta, IPenalty* p_penalty = 0, IPenalizedSolver* p_solver = 0, STK::Real threshold = STK::Arithmetic<STK::Real>::epsilon())
                      : STK::ILatentModel<STK::CArrayXX>(data),
                        y_(y),
                        beta_(beta),
                        p_penalty_(p_penalty),
                        p_solver_(p_solver),
                        n_(data->sizeRowsImpl()),
                        p_(data->sizeColsImpl()),
                        threshold_(threshold),
                        currentData_(),
                        currentBeta_(),
                        currentSet_(),
                        nbActiveVariables_(0)
          {
            initialization();
          }

      /** destructor*/
      ~PenalizedModels(){};

      //getter
      /**@return the estimated beta*/
      inline STK::CVectorX const& beta() {return beta_;}
      /** @return  pointer to the non zero beta*/
      inline STK::CVectorX const* p_currentBeta() {return &currentBeta_;}
      /** @return  the non zero beta*/
      inline STK::CVectorX const& currentBeta() {return currentBeta_;}
      /** @return  the ith non zero beta*/
      inline STK::Real const& currentBeta(int const& i) {return currentBeta_[i];}
      /** @return y the response */
      inline STK::CVectorX const& y() {return y_;}
      /** @return  a pointer to the current data*/
      inline STK::CArrayXX const* p_currentData() {return &currentData_;}
      /** @return  a pointer to the current set*/
      inline STK::Array2DVector<int> const* p_currentSet() {return &currentSet_;}
      /** @return  the current set*/
      inline STK::Array2DVector<int> const& currentSet() {return currentSet_;}
      /** @return  the ith current set*/
      inline int const& currentSet(int const& i) {return currentSet_[i];}

      //setter
      /**
       * set the solver
       * @param p_solver a pointer to a IPenalizedSolver object
       */
      inline void setSolver(IPenalizedSolver* p_solver) {p_solver_=p_solver;}
      /**
       * set the penalty of the model
       * @param p_penalty a pointer to a IPenalty object
       */
      inline void setPenalty(IPenalty* p_penalty) {p_penalty_=p_penalty;}

      /**
       * computation of the complete log-likelihood
       * @return the complete log-likelihood
       */
      STK::Real llc()
      {
        STK::Real llc = ( (y_ - (currentData_ * currentBeta_) ).square().sum() ) /2/p_penalty_->sigma2() +  p_penalty_->penaltyTerm(currentBeta_)/2;

        return llc;
      }

      /** EStep update the latent variable in penalty*/
      inline void EStep()
      {
        //STK::Real normRes=(y_- (*p_observableData_ * beta_)).norm2();//for updating sigma2_ in the penalty class

        if(p_penalty_ == 0)
        {
          std::cout<<"NO PENALTY"<<std::endl;
          this->msg_error_ = STKERROR_NO_ARG(IMultiStatModel::run(weights),parameters have not be set);
          return ;
        }
        else
          p_penalty_->update(currentBeta_);
          //p_penalty_->update(currentBeta_,normRes);
      }

      /**MStep update the estimation of beta*/
      inline void MStep()
      {
        if(p_solver_ == 0)
        {
          std::cout<<"NO SOLVER"<<std::endl;
          this->msg_error_ = STKERROR_NO_ARG(IMultiStatModel::run(weights),parameters have not be set);
          return ;
        }
        else
        {
          currentBeta_.move(p_solver_->run());
          updateCurrent();
        }

      }


      //pure virtual method from ILatentModel, no implementation here
      void SStep() {};
      void MapStep() {};

    protected:
      /**Thresholding of the new estimates : estimated coefficients < threshold_ become 0*/
      void thresholding()
      {
        for(int i = 1; i <= currentBeta_.sizeRows(); i++)
          if(std::abs(currentBeta_[i]) < threshold_)
          {
            currentBeta_[i] = 0;
            nbActiveVariables_--;
          }
      }

      /** update all the current variables*/
      void updateCurrent()
      {
        //thresholding of the new estimates
        thresholding();

        //save the new beta
        for(int i = 1; i <= currentBeta_.sizeRows(); i++)
          beta_[currentSet_[i]]=currentBeta_[i];

        if(nbActiveVariables_ != currentBeta_.sizeRows())//if TRUE, at least one variable became 0 in the previous M step
        {
          //update the currentSet
          int pos=1;
          for(int i = 1; i <= currentBeta_.sizeRows(); i++)
          {
            //erase index of the zero estimatee
            if(currentBeta_[i]==0)
            {
              currentSet_.erase(pos);
              pos--;
            }
            pos++;
          }

          //update current beta and current data
          updateCurrentData();
        }
      }

      /**Update the currentBeta_ and currentData_*/
      void updateCurrentData()
      {
        //resize currentBeta_,  currentData_ with only the active variable
        currentBeta_.resize(nbActiveVariables_);
        currentData_.resize(n_,nbActiveVariables_);
        for(int i = 1; i <= nbActiveVariables_; i++)
        {
          currentBeta_[i]=beta_[currentSet_[i]];
          currentData_.col(i)=p_observableData_->col(currentSet_[i]);
        }
      }

      /**Initialization of currentSet_, currentData_ and currentBeta_*/
      void initialization()
      {
        //thresholding of the initial solution
        for(int i = 1; i <= beta_.sizeRows(); i++)
          if(std::abs(beta_[i]) < threshold_)
            beta_[i] = 0;


        //find the index of non-zero beta
        for(int i = 1; i <= beta_.sizeRows(); i++)
        {
          if(beta_[i]!=0)
          {
            currentSet_.pushBack(1);
            currentSet_.back()=i;
          }
        }

        nbActiveVariables_ = currentSet_.sizeRows();

        //initialization of currentData_ and currentBeta_
        updateCurrentData();

      }
    private:
      ///response (length n_)
      STK::CVectorX y_;
      ///estimates of the model (length p_)
      STK::CVectorX beta_;
      ///pointer to the penalty
      IPenalty* p_penalty_;
      ///pointer to the solver for solving th M step
      IPenalizedSolver* p_solver_;
      ///size of sample
      int n_;
      ///number of covariates
      int p_;
      ///threshold under we consider a beta equal to 0
      STK::Real threshold_;
      ///dataset with only variables from currentSet
      STK::CArrayXX currentData_;
      ///beta with only non zero coefficients
      STK::CVectorX currentBeta_;
      ///index of non-zero beta coefficient
      STK::Array2DVector<int> currentSet_;
      /// number of variables in tuhe currentSet_
      int nbActiveVariables_;

  };
}

#endif /* PENALIZEDMODELS_H_ */
