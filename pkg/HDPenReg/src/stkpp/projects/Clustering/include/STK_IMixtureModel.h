/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2012  Serge Iovleff

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
 * Project:  stkpp::Clustering
 * created on: 16 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_IMixtureModel.h
 *  @brief In this file we define the main class for mixture models.
 **/


#ifndef STK_IMIXTUREMODEL_H
#define STK_IMIXTUREMODEL_H

#include "../../StatModels/include/STK_ILatentModel.h"

namespace STK
{
/** @brief Base class for Mixture model.
 *
 * In statistics, a mixture model is a probabilistic model for representing
 * the presence of sub-populations within an overall population, without
 * requiring that an observed data-set should identify the sub-population to
 * which an individual observation belongs. Formally a mixture model
 * corresponds to the mixture distribution that represents the probability
 * distribution of observations in the overall population. However, while
 * problems associated with "mixture distributions" relate to deriving the
 * properties of the overall population from those of the sub-populations,
 * "mixture models" are used to make statistical inferences about the
 * properties of the sub-populations given only observations on the pooled
 * population, without sub-population-identity information.
 *
 * Some ways of implementing mixture models involve steps that attribute
 * postulated sub-population-identities to individual observations (or weights
 * towards such sub-populations), in which case these can be regarded as types
 * unsupervised learning or clustering procedures. However not all inference
 * procedures involve such steps.
 *
 * @tparam Array can be any kind of container for the observable variables.
 * It should at least derive from ITContainer and provide an
 * access to a single row. @sa ITContainer, ICArray, IArray2DBase
 *
 * @tparam Parameters any structure encapsulating the parameters of the components.
 * @sa IMultiParameters
 **/
template <class Array, class Parameters>
class IMixtureModel : public IMixtureModelBase
{
  protected:
    using IMixtureModelBase::nbCluster_;

    /** Default constructor */
    IMixtureModel() : IMixtureModelBase(), components_() {}
    /** Constructor with data set. */
    IMixtureModel( Array const& data) : p_data_(&data), IMixtureModelBase(), components_()
    {
      this->initialize(p_data_->rows().size(), p_data_->cols().size());
      this->setRangeSamples(p_data_->cols());
    }
    /** Constructor with pointer on the data set. */
    IMixtureModel( Array const* p_data) : p_data_(p_data), IMixtureModelBase(), components_()
    {
      if (p_data_)
      {
        this->initialize(p_data_->rows().size(), p_data_->cols().size());
        this->setRangeSamples(p_data_->cols());
      }
    }
    /** copy constructor.
     *  @param model the model to copy
     **/
    IMixtureModel( IMixtureModel const& model)
                 : p_data_(model.p_data_)
                 , IMixtureModelBase(model)
                 , components_(model.components_)
    {
      for (Integer k= components_.firstIdx(); k <= components_.lastIdx(); ++k)
      { if (model.components_[k]) components_[k] = model.components_[k]->clone();}
    }
    /** destructor */
    ~IMixtureModel() {}

  public:
    /** compute the parameters. This is the default implementation
     *  of the ILatentModel::MStep abstract virtual function.
     *  This implementation assume that there is no dependencies between the
     *  parameters of the components and have to be overloaded if some parameters
     *  are shared between the components.
     **/
    virtual void mStep()
    {
      computeProportions();
      for (int k=components_.firstIdxRows(); k<= components_.lastIdxRows(); k++)
      { components_[k]->run(p_tik_->col(k));}
      this->computeLnLikelihood();
    }

    typedef IMultiStatModel<Array, Array2DVector<Real>, Parameters> IComponents;
    /** @return the array with the components */
    Array1D< IComponents* > const& components() const { return components_;}

    /** @return a constant reference on the k_th component */
    IComponents* const& p_components(int k) const { return components_[k];}
    /** @return a reference on the k_th component */
    IComponents*& p_components(int k) { return components_[k];}

    /** @return the value of the probability of the sample sample in  teh component k.
     *  @param index of the sample
     *  @param k index of the component
     **/
    virtual Real componentProbability(int i, int k)
    { return std::exp(components_[k]->computeLnLikelihood(p_data_->row(i)));}

  protected:
    /** pointer on the data set */
    Array const* p_data_;
    /** Array of the components of the mixture model */
    Array1D< IComponents* > components_;

    /** create the components */
    void createComponentsArray() { components_.resize(nbCluster_);}
    /** compute the proportions
     *  @note This method have to be overloaded if we estimate a model with
     *  equal proportions.
     **/
    virtual void computeProportions()
    {
      for (int k=this->p_prop_->firstIdx(); k<= this->p_prop_->lastIdx(); k++)
      { this->p_prop_->elt(k) = p_tik_->col(k).sum()/this->nbSample();}
    }
};

} // namespace STK

#endif /* STK_IMIXTUREMODEL_H */
