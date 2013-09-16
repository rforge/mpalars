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

/** @file STK_IMixtureModelBase.h
 *  @brief In this file we implement the abstract base class for mixture models.
 **/

#include <cmath>
#include "../include/STK_IMixtureModelBase.h"
#include "../../STatistiK/include/STK_Law_Categorical.h"
namespace STK
{
/* default constructor */
IMixtureModelBase::IMixtureModelBase() : IModelBase()
                                       , rangeSamples_(), nbCluster_(0)
                                       , p_prop_(0), p_tik_(0), p_zi_(0)
                                       , isParametersCreated_(false)
{}
/* copy constructor */
IMixtureModelBase::IMixtureModelBase( IMixtureModelBase const& model)
                                    : IModelBase(model)
                                    , rangeSamples_(model.rangeSamples_), nbCluster_(model.nbCluster_)
                    , p_prop_((isParametersCreated_ && model.p_prop_) ? model.p_prop_->clone() : model.p_prop_)
                    , p_tik_((isParametersCreated_ && model.p_tik_) ? model.p_tik_->clone() : model.p_tik_)
                    , p_zi_((isParametersCreated_ && model.p_zi_) ? model.p_zi_->clone() : model.p_zi_)
                    , isParametersCreated_(model.isParametersCreated_)
{}
/* destructor */
IMixtureModelBase::~IMixtureModelBase()
{
  if (isParametersCreated_)
  {
    if (p_prop_) delete p_prop_;
    if (p_tik_) delete p_prop_;
    if (p_zi_) delete p_prop_;
  }
}

/* initialize randomly the labels zi of the model */
void IMixtureModelBase::classInit()
{
  Law::Categorical law(*p_prop_);
  for (int i=p_zi_->firstIdxRows(); i<= p_zi_->lastIdxRows(); i++)
  { p_zi_->elt(i) = law.rand();}
  computeZik();
  mStep();
}

/* initialize randomly the posterior probabilities tik of the model */
void IMixtureModelBase::fuzziInit()
{
  RandBase generator;
  for (int i=p_tik_->firstIdxRows(); i<= p_tik_->lastIdxRows(); i++)
  {
    // create a reference on the i-th row
    Array2DPoint<Real> currentProb(p_tik_->row(i), true);
    for (Integer k= p_prop_->firstIdx(); k <= p_prop_->lastIdx(); ++k)
    {
      generator.randUnif(currentProb);
      currentProb = currentProb * (*p_prop_);
      currentProb /= currentProb.sum();
    }
  }
}

/* ceStep */
void IMixtureModelBase::ceStep()
{
  // compute tik and zi
  eStep();
  mapStep();
  // hard classification
  computeZik();
}

/* simulate zi, default implementation. */
void IMixtureModelBase::seStep()
{
  // compute tik and simulate zi
  eStep();
  for (int i=p_zi_->firstIdxRows(); i<= p_zi_->lastIdxRows(); i++)
  { p_zi_->elt(i) = Law::Categorical::rand(p_tik_->row(i));}
  // hard classification
  computeZik();
}
/* compute Tik, default implementation. */
void IMixtureModelBase::eStep()
{
  for (int i=p_tik_->firstIdxRows(); i<= p_tik_->lastIdxRows(); i++)
  {
    // compute log-likelihood of each sample for each component
    for (int k=p_tik_->firstIdxCols(); k<= p_tik_->lastIdxCols(); k++)
    { p_tik_->elt(i,k) = p_prop_->elt(k) * componentProbability(i,k);}
    // normalize
    p_tik_->row(i) /= p_tik_->row(i).sum();
  }
}
/* Compute Zi using the Map estimator, default implementation. */
void IMixtureModelBase::mapStep()
{
  for (int i=p_tik_->firstIdxRows(); i<= p_tik_->lastIdxRows(); i++)
  {
    int k;
    p_tik_->row(i).maxElt(k);
    p_zi_->elt(i) = k;
  }
}

/* compute the ln-likelihood of the mixture model. */
void IMixtureModelBase::computeLnLikelihood()
{
  Real sum = 0.;
  for (int i = rangeSamples_.firstIdx(); i<= rangeSamples_.lastIdx(); ++i)
  {
    for (int k=p_prop_->firstIdx(); k<= p_prop_->lastIdx(); k++)
    { sum += p_prop_->elt(k) * componentProbability(i,k);}
  }
  setLnLikelihood(std::log(sum));
}

/* set the parameters of the  mixture model.
 *  @param p_prop pointer on the proportion of the mixture model
 *  @param p_tik pointer on the posterior probabilities
 *  @param p_zi pointer on the class labels
 * */
void IMixtureModelBase::setMixtureParameters( Array2DPoint<Real>* p_prop
                                            , Array2D<Real>* p_tik
                                            , Array2DVector<int>* p_zi)
{
  if (isParametersCreated_)
  {
    if (p_prop_) delete p_prop_;  p_prop_ = 0;
    if (p_tik_) delete p_prop_;  p_tik_ = 0;
    if (p_zi_) delete p_prop_;  p_zi_ = 0;
  }
  p_prop_ = p_prop;
  p_tik_ = p_tik;
  p_zi_ = p_zi;
  isParametersCreated_ = false;
}

/* Create the parameters of the  mixture model. */
void IMixtureModelBase::createMixtureParameters()
{
  createProp();
  createTik();
  createZi();
  isParametersCreated_ = true;
}

/* create the proportions */
void IMixtureModelBase::createProp()
{
  if (p_prop_) { p_prop_->resize(nbCluster_);}
  else         { p_prop_ = new Array2DPoint<Real>(nbCluster_);}
  p_prop_->setValue(1./nbCluster_);
}
/* create the tik probabilities */
void IMixtureModelBase::createTik()
{
  if (p_tik_) { p_tik_->resize(rangeSamples_, nbCluster_);}
  else        { p_tik_ = new Array2D<Real>(rangeSamples_, nbCluster_);}
  p_tik_->setValue(1./nbCluster_);
}
/* create the zi labels */
void IMixtureModelBase::createZi()
{
  if (p_zi_) { p_zi_->resize(rangeSamples_);}
  else       { p_zi_ = new Array2DVector<int>(rangeSamples_);}
  p_zi_->setValue(STKBASEARRAYS);
}

/* replace tik by zik, the indicator variable of the zi */
void IMixtureModelBase::computeZik()
{
  (*p_tik_) = 0.;
  for (int i=p_tik_->firstIdxRows(); i<= p_tik_->lastIdxRows(); i++)
  { p_tik_->elt(i, p_zi_->elt(i)) = 1.;}
}

} // namespace SDTK
