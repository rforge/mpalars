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
 * created on: 23 mai 2013
 * Author:   Quentin Grimonprez
 **/

/** @file PenalizedModels.h
 *  @brief In this file, we find the definition of the PenalizedModels class.
 **/


#ifndef PENALIZEDMODELS_H_
#define PENALIZEDMODELS_H_

#include "../../stkpp/projects/StatModels/include/STK_IModelBase.h"
#include "../../stkpp/include/STKpp.h"
#include "IPenalizedSolver.h"
#include "IPenalty.h"

namespace HD
{
  /** @ingroup lassoModels
   *  @brief The class PenalizedModels inherits from the template class @c ILatentModel.
   *  It contains implementations of Estep and Mstep for run an EM algorithm.
   */

  class PenalizedModelsBase : public STK::IModelBase
  {
    public:
      /**
       * Constructor
       * @param p_data each row contains values of covariates for an individuals
       * @param y response of each individuals
       * @param beta Initialization of beta_ (least-Square solution)
       * @param p_penalty pointer on the penalty of the model
       * @param p_solver pointer on the solver to use in the Mstep
       * @param threshold tolerance for thresholding the solution
       */
      PenalizedModelsBase(STK::CArrayXX const* p_data, STK::CVectorX const* p_y, STK::CVectorX const& beta, IPenalizedSolver* p_solver = 0)
                      : STK::IModelBase(p_data->sizeRows(),p_data->sizeCols()),
                        p_data_(p_data),
                        p_y_(p_y),
                        beta_(beta)
          {
          }

      PenalizedModelsBase(STK::CArrayXX const* p_data, STK::CVectorX const* p_y)
                      : STK::IModelBase(p_data->sizeRows(),p_data->sizeCols()),
                        p_data_(p_data),
                        p_y_(p_y),
                        beta_(p_y->sizeRows())
          {
          }

      PenalizedModelsBase()
                      : STK::IModelBase(),p_data_(0),p_y_(0),beta_()
          {}

      /** destructor*/
      virtual ~PenalizedModelsBase(){};

      /** copy constructor */
      //PenalizedModels

      /** clone */
      //PenalizedModels* clone() const { return new();};

      //setter
      inline void setP_y(STK::CVectorX const* p_y) {p_y_ = p_y;}
      inline void setP_data(STK::CArrayXX const* p_data)
      {
        p_data_ = p_data;
        setNbSample(p_data_->sizeRows());
        setNbVar(p_data_->sizeCols());
      }
      //getter
      /**@return a pointer to the data*/
      inline STK::CArrayXX const* p_data() const {return p_data_;}
      /**@return the estimated beta*/
      inline STK::CVectorX const& beta() const {return beta_;}
      /**@return the estimated beta*/
      inline STK::Real beta(int i) const {return beta_[i];}
      /** @return a pointer to the response */
      inline STK::CVectorX const* p_y() const {return p_y_;}

      /** EStep update the latent variable in penalty*/
      virtual void eStep(){}

      /**MStep update the estimation of beta*/
      virtual void mStep(bool burn = true){}


    protected:
      ///pointer to the data
      STK::CArrayXX const* p_data_;
      ///response (length n_)
      STK::CVectorX const* p_y_;
      ///estimates of the model (length p_)
      STK::CVectorX beta_;

  };


  template <class Model> struct ModelTraits;

  template<class Model>
  class PenalizedModels : public PenalizedModelsBase
  {
    public:

    typedef typename ModelTraits<Model>::Solver Solver;
    typedef typename ModelTraits<Model>::Penalty Penalty;

      /**
       * Constructor
       * @param p_data each row contains values of covariates for an individuals
       * @param y response of each individuals
       * @param beta Initialization of beta_ (least-Square solution)
       * @param p_penalty pointer on the penalty of the model
       * @param p_solver pointer on the solver to use in the Mstep
       * @param threshold tolerance for thresholding the solution
       */
      PenalizedModels(STK::CArrayXX const* p_data, STK::CVectorX const* p_y, STK::CVectorX const& beta, Solver* p_solver = 0)
                      : PenalizedModelsBase(p_data,p_y,beta),
                        p_solver_(p_solver)
          {
          }

      PenalizedModels(STK::CArrayXX const* p_data, STK::CVectorX const* p_y)
                      : PenalizedModelsBase(p_data,p_y),
                        p_solver_(0)
          {
          }

      PenalizedModels()
                      : PenalizedModelsBase(),p_solver_(0)
          {}

      /** destructor*/
       virtual ~PenalizedModels()
       {
         if(p_penalty_) delete p_penalty_;
         if(p_solver_) delete p_solver_;
       };


      //getter
      /** @return  the non zero beta*/
      inline STK::CVectorX const& currentBeta() const {return p_solver_->currentBeta();}
      /** @return  the non zero beta*/
      inline STK::CVectorX const* p_currentBeta() const {return p_solver_->p_currentBeta();}
      /** @return  data*/
      inline STK::CArrayXX const& currentData() const {return p_solver_->currentData();}
      /** @return  data*/
      inline STK::CArrayXX const* p_currentData() const {return p_solver_->p_currentData();}
      /**@return currentSet_*/
      inline STK::Array2DVector<int> const& currentSet() const {return p_solver_->currentSet();};
      /**@return p_solver_*/
      inline Solver const* p_solver() const {return p_solver_;}
      /**@return p_penalty_*/
      inline Penalty const* p_penalty() const {return p_solver_->p_penalty();}

      //setter
      /**
       * set the solver
       * @param p_solver a pointer to a IPenalizedSolver object
       */
      inline void setSolver(Solver* p_solver) {p_solver_=p_solver;}



      /** EStep update the latent variable in penalty*/
      inline void eStep()
      {
        //STK::Real normRes=(y_- (*p_observableData_ * beta_)).norm2();//for updating sigma2_ in the penalty class

        if(p_solver_ == 0)
          throw STK::invalid_argument(STK::String("p_solver_ has not be set"));
        else
          p_solver_ ->update();

      }

      /**MStep update the estimation of beta*/
      inline void mStep(bool burn = true)
      {
        if(p_solver_ == 0)
          throw STK::invalid_argument(STK::String("p_solver_ has not be set"));
        else
        {
          setLnLikelihood(p_solver_->run(burn));
          beta_=p_solver_->beta();
        }
      }

    protected:
      ///pointer to the solver for solving the M step
      Solver* p_solver_;
      Penalty* p_penalty_;


  };
}

#endif /* PENALIZEDMODELS_H_ */
