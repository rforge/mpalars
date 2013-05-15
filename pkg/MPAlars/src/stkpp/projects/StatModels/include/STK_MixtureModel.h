
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

/** @file STK_MixtureModel.h
 *  @brief In this file we define the main class fort mixture models.
 **/


#ifndef STK_MIXTUREMODEL_H
#define STK_MIXTUREMODEL_H

#include "STK_LatentModel.h"

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
 **/
template <class ObservableData>
class MixtureModel : public LatentModel<ObservableData, Array2D< Real > >
{
  public:
    typedef LatentModel<ObservableData, Array2D<Real> > Base;

  protected:
    /** Constructor with data set. */
    MixtureModel(ObservableData const& data) : Base(data) {}
      /** Constructor with data set. */
    MixtureModel(ObservableData const* p_data) : Base(p_data) {}
    /** */
    Array2D<Real> const* p_tik_() const { return this->p_latentData();}

  public:
    /** destructor */
    virtual ~MixtureModel() { ;}

};

} // namespace STK

#endif /* STK_MIXTUREMODEL_H */
