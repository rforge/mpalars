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
 * Project:  stkpp::
 * created on: 16 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_IMixtureModel.h
 *  @brief In this file we define the main class fort mixture models.
 **/


#ifndef STK_IMIXTUREMODEL_H
#define STK_IMIXTUREMODEL_H

#include "STK_ILatentModel.h"

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
 * In this interface we assume there is an underline generative model that will
 * be estimated using either an EM, SEM or CEM algorithm.
 *
 * @tparam Array can be any kind of container for the observable variables.
 * It should at least derive from ITContainer and provide an
 * access to a single row. @sa ITContainer, ICArray, IArray2DBase
 *
 * @tparam Parameters any structure encapsulating the parameters of the components.
 * @sa IMultiParameters
 *
 **/
template <class Array, class Parameters>
class IMixtureModel : public ILatentModel<Array >
{
  protected:
    typedef ILatentModel<Array > Base;
    using Base::p_data_;

    /** Default constructor */
    IMixtureModel() : Base(), p_prop_(0), p_tik_(0), p_zi_(0), components_() {}
    /** Constructor with data set. */
    IMixtureModel(Array const& data) : Base(data), p_prop_(0), p_tik_(0), p_zi_(0), components_()
    {}
      /** Constructor with data set. */
    IMixtureModel(Array const* p_data) : Base(p_data), p_prop_(0), p_tik_(0), p_zi_(0), components_()
    {}
    /** destructor */
    ~IMixtureModel() {}

  public:
    typedef IMultiStatModel<Array, Array2DVector<Real>, Parameters> IComponents;
    /** @return the number of cluster. */
    int nbCluster() const { return nbCluster_;}
    /** @return the proportions of each mixtures */
    Array2DPoint<Real> const* p_prop() const { return p_prop_;}
    /** @return a constant pointer on the tik probabilities */
    Array2D<Real> const* p_tik() const { return p_tik_;}
    /** @return a constant pointer on the zi labels */
    Array2D<int> const* p_zi() const { return p_zi_;}
    /** @return the array with the components */
    Array1D< IComponents* > const& components() const { return components_;}

    /** @return a constant reference on the k_th component */
    IComponents* const& p_component(int k) const { return components_[k];}
    /** @return a reference on the k_th component */
    IComponents*& p_component(int k) { return components_[k];}

    /** set the proportions */
    void setNbCluster(int nbCluster) { nbCluster_ = nbCluster;}
    /** set the proportions */
    void setProp(Array2D<Real> const& prop) { p_prop_ = &prop;}
    /** set the tik probabilities */
    void setTik(Array2D<Real> const& tik) { p_tik_ = &tik;}
    /** set the zi labels */
    void setZi(Array2D<int> const& zi) { p_zi_ = &zi;}

    /** create the proportions */
    void createProp()
    {
      if (p_prop_) { p_prop_->resize(nbCluster);}
      else         { p_prop_ = new Array2DPoint<Real>(nbCluster);}
      p_prop_->setValue(1./nbCluster_);
    }
    /** create the tik probabilities */
    void createTik()
    {
      if (p_tik_) { p_tik_->resize(p_data_->cols(), nbCluster);}
      else        { p_tik_ = new Array2D<Real>(p_data_->cols(), nbCluster);}
      p_tik_->setValue(1./nbCluster_);
   }
    /** create the zi labels */
    void createZi()
    {
      if (p_zi_) { p_zi_->resize(p_data_->cols());}
      else       { p_zi_ = new Array2DVector<Real>(p_data_->cols());}
      MapStep();
    }
    /** create the components */
    void createComponents()
    {
      components_.resize(nbCluster_);
      for (int k=components_.firstIdxRows(); k<= components_.lastIdxRows(); k++)
      {
        if (!components_[k]) { components_[k] = new IComponents();}
        components_[k]->setData(p_data_);
      }
    }

    /** Compute Z using the Map estimator. This is the concrete implementation
     *  of the ILatentModel::MapStep abstract virtual function.
     **/
    virtual void MapStep()
    {
      if (!p_tik_)
      {STKRUNTIME_ERROR_NO_ARG(IMixtureModel::MapStep,Tik is not set nor created);}
      if (!p_zi_)
      {STKRUNTIME_ERROR_NO_ARG(IMixtureModel::MapStep,Zi is not set nor created);}

      for (int i=p_tik_->firstIdxRows(); i<= p_tik_->lastIdxRows(); i++)
      {
        int k;
        p_tik_->row(i).maxElt(k);
        p_zi_->elt(i) = k;
      }
    }
    /** Compute Z using the Map estimator. */
    inline virtual void CStep() { MapStep();}
    /** compute Tik. This is the concrete implementation
     *  of the ILatentModel::EStep abstract virtual function.
     **/
    virtual void EStep()
    {
      if (!p_data_)
      {STKRUNTIME_ERROR_NO_ARG(IMixtureModel::EStep,data is not set);}
      if (!p_tik_)
      {STKRUNTIME_ERROR_NO_ARG(IMixtureModel::EStep,Tik is not set nor created);}

      for (int i=p_tik_->firstIdxRows(); i<= p_tik_->lastIdxRows(); i++)
      {
        // compute log-likelihood of each sample for each component
        for (int k=p_tik_->firstIdxCols(); k<= p_tik_->lastIdxCols(); k++)
        { p_tik_->elt(i,k) = std::exp(components_[k]->computelnLikelihood(p_data_->row(i)));}
        // normalize
        p_tik_->row(i) /= p_tik_->row(i).sum();
      }
    }
    /** compute the parameters. This is the concrete implementation
     *  of the ILatentModel::MStep abstract virtual function.
     *  This implementation assume that there is no dependencies between the
     *  parameters of the components.
     **/
    virtual void MStep()
    {
      if (!p_data_)
      {STKRUNTIME_ERROR_NO_ARG(IMixtureModel::EStep,data is not set);}
      if (!p_tik_)
      {STKRUNTIME_ERROR_NO_ARG(IMixtureModel::EStep,Tik is not set nor created);}

      for (int k=components_.firstIdxRows(); k<= components_.lastIdxRows(); k++)
      { components_[k]->run(p_tik_->col(k));}
    }

  protected:
    /** number of cluster. */
    int nbCluster_;
    /** The proportions of each mixtures */
    Array2DPoint<Real>* p_prop_;
    /** The tik probabilities */
    Array2D<Real>* p_tik_;
    /** The zik class label */
    Array2DVector<int>* p_zi_;
    /** Array of the components of the mixture model */
    Array1D< IComponents* > components_;
};

} // namespace STK

#endif /* STK_IMIXTUREMODEL_H */
