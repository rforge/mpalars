/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2011  Serge Iovleff

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

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  stkpp::Model
 * created on: 22 juil. 2011
 * Purpose: define the Interface base class LatentModel (Statistical Model).
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_LatentModel.h
 *  @brief In this file we define the Interface base class LatentModel.
 **/

#ifndef STK_LATENTMODEL_H
#define STK_LATENTMODEL_H

#include <cmath>

#include "STK_StatisticalModel.h"
#include "../../STatistiK/include/STK_Law_IMultiLaw.h"

namespace STK
{

/** @ingroup StatModels
 *  @brief Base class for latent models.
 *  In statistics, latent variables (as opposed to observable variables)
 *  are variables that are not directly observed but are rather inferred
 *  (through a mathematical model) from other variables that are observed
 *  (directly measured). Mathematical models that aim to explain observed
 *  variables in terms of latent variables are called latent variable models.
 *
 *  Generally it is assumed that:
 *  -the responses on the observable variables are the result of an
 *   individual's position on the latent variable(s) and,
 * - that the observable variables have nothing in common after controlling for
 *   the latent variable (local independence).
 *
 *  Sometimes latent variables correspond to aspects of physical reality,
 *  which could in principle be measured, but may not be for practical
 *  reasons. In this situation, the term hidden variables is commonly used
 *  (reflecting the fact that the variables are "really there", but hidden).
 *  Other times, latent variables correspond to abstract concepts, like
 *  categories, behavioral or mental states, or data structures.
 *
 *  This class propose some virtual function that have to be re-implemented by
 *  derived classes.
 *
 *  @tparam ArrayObservable can be any kind of container for the observable variables.
 *  It should at least derive from IContainer2D and provide an
 *  access to a single row. @sa IContainer2D
 *  @tparam ArrayLatent can be any kind of container for the latent variables.
 *  It should at least derive from IContainer2D and provide an
 *  access to a single row. @sa IContainer2D
 **/
template <class ArrayObservable, class ArrayLatent>
class LatentModel : public IModelBase, public IRunnerBase
{
  protected:
    /** Constructor with data set
     *  @param data the observable data set
     **/
    LatentModel(ArrayObservable const& data) : IModelBase(data.sizeRows(), data.sizeCols())
                                             , IRunnerBase()
                                             , p_latentData_(0)
    {}
    /** Constructor with a pointer on the data set
     *  @param p_data a pointer on the constant data set
     **/
    LatentModel(ArrayObservable const* p_data) : p_observableData_(p_data)
                                               , IRunnerBase()
                                               , p_latentData_(0)
    {
      if (p_observableData_)
      { this->initialize(p_observableData_->sizeRows(), p_observableData_->sizeCols());}
    }
  public:
    /** destructor */
    virtual ~LatentModel() {}
    /** @return a pointer on the latent data set*/
    ArrayObservable const* p_data() const { return p_observableData_;}
    /** Type of the row container (the sample) */
    typedef typename ArrayObservable::Row Row;
    /** Expectation step of the EM algorithm. */
    virtual void EStep() =0;
    /** Maximization step of the EM algorithm. */
    virtual void MStep() =0;
    /** simulation step of the SEM algorithm. */
    virtual void SStep() =0;
    /** Posterior maximization of the latent variable. In a classification model,
     *  this is the  Classification step of the CEMalgorithm.
     **/
    virtual void MapStep() =0;

  protected:
    /** By default the virtual method @c run inherited from IRunnerBase does
     *  nothing */
    virtual bool run() { return true;}
    /** A pointer on the original data set */
    ArrayObservable const* p_observableData_;
    /** The latent data set */
    ArrayLatent* p_latentData_;
};

} // namespace STK

#endif /* STK_IMODEL_H */
