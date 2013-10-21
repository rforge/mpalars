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
 * created on: 18 sept. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file FusedLassoSolver.cpp
 *  @brief In this file, implementation of the methods of the @c FusedLassoSolver class. .
 **/

#include "FusedLassoSolver.h"

namespace HD
{
  /*default constructor*/
  FusedLassoSolver::FusedLassoSolver()
  : IPenalizedSolver()
  , p_solver_(0)
  , currentXty_()
  , segment_()
  , n_(0)
  , nbActiveVariables_(0)
  , eps_(STK::Arithmetic<STK::Real>::epsilon())
  , p_penalty_(0)
  {
  }

  /*
   * Constructor
   * @param p_data pointer to the full data
   * @param beta initial solution of the problem
   * @param p_y pointer to the response associated with the data
   * @param burn burn-in period before regrouping variables in segments
   * @param p_solver pointer to the solver
   * @param p_penalty pointer to the Fused lasso penalty
   * @param eps tolerance for zero
   */
  FusedLassoSolver::FusedLassoSolver(STK::CArrayXX const* p_data, STK::CVectorX const& beta, STK::CVectorX const* p_y,
                                     STK::CG<FusedLassoMultiplicator,STK::CVectorX,InitFunctor>* p_solver, FusedLassoPenalty* p_penalty,
                                     STK::Real eps)
                                     : IPenalizedSolver(beta, p_data, p_y)
                                     , p_solver_(p_solver)
                                     , currentXty_()
                                     , segment_(beta.sizeRows())
                                     , n_(p_data->sizeRows())
                                     , nbActiveVariables_(beta.sizeRows())
                                     , eps_(eps)
                                     , p_penalty_(p_penalty)
  {
    //initialization of currentXty
    currentXty_ = p_data->transpose() * (*p_y);

    p_solver_->setB(currentXty_);

    for(int i = 1; i <= beta.sizeRows(); i++)
    {
      //initialization of the segment : segment of 1 point
      segment_[i-1] = STK::Range(i,1);
    }
  }


  /*update the penalty (E step)*/
  void FusedLassoSolver::update()
  {
    if(p_penalty_ == 0)
      throw STK::invalid_argument(STK::String("p_penalty_ has not be set"));

    p_penalty_->update(currentBeta_);
  }

  /*initialize the container of the class*/
  void FusedLassoSolver::initializeSolver()
  {
    //check the existence pointers to the data and the response
    if(p_data_ == 0)
      throw STK::invalid_argument(STK::String("p_data_ has not be set"));
    if(p_y_ == 0)
      throw STK::invalid_argument(STK::String("p_y_ has not be set"));

    currentData_ = *p_data_;
    currentBeta_ = beta_;
    nbActiveVariables_ = p_data_->sizeCols();
    n_ = p_data_->sizeRows();
//    currentXty_.resize(nbActiveVariables_);
    currentXty_ = currentData_.transpose() * (*p_y_);

    p_solver_->setB(currentXty_);
    currentSet_.resize(nbActiveVariables_);

    segment_.resize(nbActiveVariables_);

    for(int i = 1; i <= nbActiveVariables_; i++)
    {
      //initialization of the segment : segment of 1 point
      segment_[i-1] = STK::Range(i,1);
      //initialization of the set
      currentSet_[i] =  i;
    }

  }

  /*Computation of the completed log-likelihood
   * @return the current completed loglikelihood
   */
  STK::Real FusedLassoSolver::computeLlc()
  {
    STK::Real llc, temp (0);
    //when regrouping variables in segments, the part : beta^(k+1) ^2/|beta^(k)| is missing for all variable in the segment
    for(int i = 0; i < (int) segment_.size(); i++)
      temp += (currentBeta_[i+1] * currentBeta_[i+1]) * (segment_[i].size()-1) / (std::abs(beta_[segment_[i].firstIdx()]) + p_penalty_->eps());
    temp *= p_penalty_->lambda1();

    llc = - ( ( (*p_y_ - (currentData_ * currentBeta_) ).square().sum() )/p_penalty_->sigma2() +  p_penalty_->penaltyTerm(currentBeta_) + temp)/2;

    return llc;
  }

  /*
   * Solve the M-step with a conjugate gradient
   * @return new estimate of beta
   */
  STK::Real FusedLassoSolver::run(bool const& burn)
  {
//    p_solver_->setB(currentXty_);

    //run the conjugate gradient
    p_solver_->run();

    //get the solution
    currentBeta_ = p_solver_->x();

    //compute llc
    STK::Real llc = computeLlc();

    //reduction of the data with segments
    if(burn)
    {
      //update segment_ and currentSet_ and currentBeta_
      bool changement = false;
      changement = updateCurrent();

      //update the full beta_
      updateBeta();

      //if true, we have to change the currentData matrix
      if(changement)
      {
        updateCurrentData();

        //update currentXty_ because currentData_ change
//        currentXty_.resize(currentData_.sizeCols());
        currentXty_ = currentData_.transpose() * (*p_y_);
      }
    }
    else
      beta_=currentBeta_;

    return llc;
  }

  /*
   * update currentSet_, segment_ and currentBeta_
   * @return a boolean, if true there is a changement in the currentSet
   */
  bool FusedLassoSolver::updateCurrent()
  {
    int betaSize = currentBeta_.sizeRows();
    STK::Array2DVector<STK::Real> betaTemp (betaSize);
    betaTemp=currentBeta_;
    bool changement = false;

    for(int i = betaSize; i > 1; i--)
    {
      //fusion of 2 segments
      if( std::abs(currentBeta_[i]-currentBeta_[i-1]) <= eps_ )
      {
        //erase the i-th beta (same as (i-1)-th)
        betaTemp.erase(i,1);

        /*update of segment_*/
        //increase the last index of the segment i-1 //vector index shift by 1
        segment_[i-2] = segment_[i-2].incLast(segment_[i-1].lastIdx()-segment_[i-1].firstIdx()+1);

        //delete the segment i
        segment_.erase(segment_.begin()+i-1,segment_.begin()+i);

        /*update of indBetaSegment*/
        //the segment i merged with segment i-1, segment i+1, i+2,... become segment i, i+1, ...
        for( int j = i; j <= currentSet_.sizeRows(); j++)
          currentSet_[j]--;

        //a variable is removed from the active variable
        nbActiveVariables_--;
        changement = true;
      }
    }

    //put the changes into currentBeta
    if(changement)
    {
//      currentBeta_.resize(nbActiveVariables_);
      currentBeta_ = betaTemp;
    }

    for(int i = 1; i <= currentBeta_.sizeRows(); i++)
    {
      if( std::abs(currentBeta_[i]) < eps_ )
        currentBeta_[i] = 0;
    }

    return changement;
  }

  /*
   * update currentData_
   */
  void FusedLassoSolver::updateCurrentData()
  {
    //resize currentData_ with only the active variable
    currentData_.resize(STK::Range(1,n_),STK::Range(1,nbActiveVariables_));//currentData_.resize(n_,nbActiveVariables_);
    currentData_.zeros();

    for(int i = 0; i < nbActiveVariables_; i++)
    {
      //the covariate associates to an index is the sum of all covariates of the segment
      for(int j = segment_[i].firstIdx(); j<= segment_[i].lastIdx(); j++)
        currentData_.col(i+1) += p_data_->col(j);
    }
  }

  /*
   * update beta_
   */
  void FusedLassoSolver::updateBeta()
  {
    //complete beta with the updated values
    for(int i = 0; i < nbActiveVariables_; i++)
    {
      for(int j = segment_[i].firstIdx(); j <= segment_[i].lastIdx(); j++)
        beta_[j] = currentBeta_[i+1];
    }
  }

}

