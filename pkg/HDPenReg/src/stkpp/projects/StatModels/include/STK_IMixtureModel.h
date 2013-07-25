
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


#ifndef STK_IMixtureModel_H
#define STK_IMixtureModel_H

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
 * In this class we assume there is a underline generative model that will be
 * estimated using either an EM, SEM or CEM algorithm.
 **/
template <class Array, class Components>
class IMixtureModel : public ILatentModel<Array >
{
  public:
    typedef ILatentModel<Array > Base;
    /** destructor */
    virtual ~IMixtureModel() : Base(), p_tik_(0), p_zi_(0), components_() {}
    /** Constructor with data set. */
    IMixtureModel(Array const& data) : Base(data) {}
      /** Constructor with data set. */
    IMixtureModel(Array const* p_data) : Base(p_data) {}
    /** @return a pointer on the tik probabilities */
    Array2D<Real> const* p_tik() const { return p_tik_;}
    /** @return a constant pointer on the zi labels */
    Array2D<int> const* p_zi() const { return p_zi_;}
    /** set the t_ik probabilities */
    void set_tik(Array2D<Real> const& t_ik) { p_tik_ = &t_ik;}
    /** set the zi labels */
    void set_zi(Array2D<int> const& zi) { p_zi_ = &zi;}

  protected:
    /** The proportions of each mixtures */
    Array2DPoint<Real>* p_prop_;
    /** The tik probabilities */
    Array2D<Real>* p_tik_;
    /** The zik class label */
    Array2DVector<int>* p_zi_;
    /** Array of the components of the models */
    Array2DPoint<Components*> components_;
    /** comput zi using a Map estimator using the tik */
    void computeZ()
    {
      for (int i=p_tik_->firstIdxRows(); i<= p_tik_->lastIdxRows(); i++)
      {
        int k;
        p_tik_->row(i).maxElt(k);
        p_zi_->elt(i) = k;
      }
    }
};

} // namespace STK

#endif /* STK_IMixtureModel_H */
