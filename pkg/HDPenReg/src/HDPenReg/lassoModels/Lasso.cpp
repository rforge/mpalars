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
 * created on: 17 déc. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file Lasso.cpp
 *  @brief In this file .
 **/

#include "Lasso.h"


namespace HD
{
  /* default constructor*/
  Lasso::Lasso()
  : PenalizedModels<Lasso>()
  {
    STK::CVectorX Xty(nbVariable());

    // creation lasso penalty
    LassoPenalty*  p_penalty = new LassoPenalty();
    p_penalty_ = p_penalty;

    //creation functor for CG
    LassoMultiplicator mult(p_currentData(),p_penalty_->p_invPenalty(), p_penalty_->p_sigma2());
    mult_ = mult;

    //create CG
    STK::CG<LassoMultiplicator,STK::CVectorX, InitFunctor>* p_gcsolver = new STK::CG<LassoMultiplicator,STK::CVectorX,InitFunctor>(mult_, Xty);
    p_gcsolver_ = p_gcsolver;

    //create solver for lasso
    LassoSolver* p_lassosolver = new LassoSolver();
    p_lassosolver->setPenalty(p_penalty_);
    p_lassosolver->setSolver(p_gcsolver_);

    //add solver to lasso
    p_solver_ = p_lassosolver;

    //add pointer to currentData to the multiplicator functor
    mult_.p_data_= p_currentData();
  }

  /*
   * Constructor
   * @param p_data pointer to the data
   * @param p_y pointer to the response
   * @param lambda value of parameter associated to the l1 penalty
   * @param threshold threshold for setting coefficient to 0
   * @param epsCG epsilon for CG convergence
   */
  Lasso::Lasso(STK::CArrayXX const* p_data, STK::CVectorX const* p_y, STK::Real lambda, STK::Real threshold, STK::Real epsCG)
  : PenalizedModels<Lasso>(p_data, p_y)
  {
    STK::CVectorX Xty(nbVariable());
    Xty = p_data_->transpose() * *p_y_;

    STK::CVectorX beta0(Xty.sizeRows());
    for(int i = 1; i <= Xty.sizeRows();i++)
      beta0[i] = Xty[i]/((*p_data).col(i).norm2());

    beta_ = beta0;

    // creation lasso penalty
    LassoPenalty*  p_penalty = new LassoPenalty(lambda);
    p_penalty_ = p_penalty;
    //creation functor for CG
    LassoMultiplicator mult(p_currentData(),p_penalty_->p_invPenalty(), p_penalty_->p_sigma2());
    mult_ = mult;

    //create CG
    STK::CG<LassoMultiplicator,STK::CVectorX, InitFunctor>* p_gcsolver = new STK::CG<LassoMultiplicator,STK::CVectorX,InitFunctor>(mult_, Xty, 0, epsCG);
    p_gcsolver_ = p_gcsolver;

    //create solver for lasso
    LassoSolver* p_lassosolver = new LassoSolver(p_data_, beta_, p_y_, threshold, p_gcsolver_, p_penalty_);

    //add the penalty
    p_lassosolver->setPenalty(p_penalty_);

    //add solver to lasso
    p_solver_ = p_lassosolver;

//        InitFunctor init(p_solver_->p_currentBeta());
//        init_ = init;
//        (p_solver_->p_solver())->setInitFunctor(&init_);

    mult_.p_data_= p_currentData();
  }

  /* destructor*/
  Lasso::~Lasso()
  {
    if(p_gcsolver_) delete p_gcsolver_;
  }

  /*initialize the containers of all subclasses*/
  void Lasso::initializeModel()
  {
    //compte intitial beta
    STK::CVectorX Xty(nbVariable());
    Xty = p_data_->transpose() * *p_y_;

    STK::CVectorX beta0(Xty.sizeRows());
    for(int i = 1; i <= Xty.sizeRows();i++)
      beta0[i] = Xty[i]/((*p_data_).col(i).norm2());

    beta_ = beta0;
    //set the parameter of the sover
    p_solver_->setData(p_data_);
    p_solver_->setY(p_y_);
    p_solver_->setBeta(beta_);
    //initialize the solver
    p_solver_->initializeSolver();
  }

  /*initialization of the class for with a new beta0
   * @param beta initial start for beta
   * */
  void Lasso::initializeBeta(STK::CVectorX const& beta)
  {
    p_solver_->setBeta(beta_);
    p_solver_->initializeSolver();
  }
}


