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
#include "../../Arrays/include/STK_Array2D.h"
#ifdef STK_MIXTURE_DEBUG
#include "../../Arrays/include/STK_Display.h"
#endif
#include "../include/STK_IMixtureModelBase.h"
#include "../../STatistiK/include/STK_Law_Categorical.h"

namespace STK
{
/* default constructor */
IMixtureModelBase::IMixtureModelBase(int nbCluster) : IModelBase()
                                       , nbCluster_(nbCluster)
                                       , p_prop_(0), p_tik_(0), p_zi_(0)
                                       , isParametersCreated_(false)
                                       , state_(Clust::modelCreated_)
{}
/* copy constructor */
IMixtureModelBase::IMixtureModelBase( IMixtureModelBase const& model)
                                    : IModelBase(model)
                                    , nbCluster_(model.nbCluster_)
                    , p_prop_((model.isParametersCreated_ && model.p_prop_) ? model.p_prop_->clone() : model.p_prop_)
                    , p_tik_((model.isParametersCreated_ && model.p_tik_) ? model.p_tik_->clone() : model.p_tik_)
                    , p_zi_((model.isParametersCreated_ && model.p_zi_) ? model.p_zi_->clone() : model.p_zi_)
                    , isParametersCreated_(model.isParametersCreated_)
                    , state_(Clust::modelCreated_)
{}
/* destructor */
IMixtureModelBase::~IMixtureModelBase()
{ deleteMixtureParameters();}

/* This function can be overloaded in derived class for initialization of
 *  the mixture parameters.
 **/
void IMixtureModelBase::initializeModel()
{
  // create mixture parameters
  if (!isParametersCreated_ && !(p_prop_ && p_tik_ && p_zi_)) { createMixtureParameters();}
  this->setNbFreeParameters(computeNbFreeParameters());
  state_ = Clust::modelInitialized_;
}

/* initialize randomly the labels zi of the model */
void IMixtureModelBase::randomClassInit()
{
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Entering IMixtureModelBase::randomClassInit()\n");
#endif
  if (isParametersCreated_) *p_prop_ = 1./Real(nbCluster_);
  Law::Categorical law(*p_prop_);
  for (int i = p_zi_->firstIdx(); i<= p_zi_->lastIdx(); ++i)
  { p_zi_->elt(i) = law.rand();}
  cStep();
  initializeStep();
  eStep();
}

/* initialize randomly the posterior probabilities tik of the model */
void IMixtureModelBase::randomFuzzyInit()
{
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Entering IMixtureModelBase::randomFuzzyInit()\n");
#endif
  if (isParametersCreated_) *p_prop_ = 1./Real(nbCluster_);
  RandBase generator;
  for (int i = p_tik_->firstIdxRows(); i<= p_tik_->lastIdxRows(); ++i)
  {
    // create a reference on the i-th row
    Array2DPoint<Real> tikRowi(p_tik_->row(i), true);
    generator.randUnif(tikRowi);
    tikRowi = tikRowi * (*p_prop_);
    tikRowi /= tikRowi.sum();
  }
  initializeStep();
  eStep();
}

/* cStep */
void IMixtureModelBase::cStep()
{
  (*p_tik_) = 0.;
  for (int i=p_tik_->firstIdxRows(); i<= p_tik_->lastIdxRows(); i++)
  { p_tik_->elt(i, p_zi_->elt(i)) = 1.;}
}

/* simulate zi  */
void IMixtureModelBase::sStep()
{
  // simulate zi
  for (int i = p_zi_->firstIdx(); i<= p_zi_->lastIdx(); ++i)
  { p_zi_->elt(i) = Law::Categorical::rand(p_tik_->row(i));}
  // hard classification
  cStep();
}
/* compute Tik, default implementation. */
void IMixtureModelBase::eStep()
{
  Real sum = 0.;
  for (int i = p_tik_->firstIdxRows(); i<= p_tik_->lastIdxRows(); ++i)
  {
    Array2DPoint<Real> lnComp(p_tik_->cols());
    for (int k=p_tik_->firstIdxCols(); k<= p_tik_->lastIdxCols(); k++)
    { lnComp[k] = lnComponentProbability(i,k);}
    int kmax;
    Real max = lnComp.maxElt(kmax);
    p_zi_->elt(i) = kmax;
    // compute sum_k pk exp{lnCom_k - lnComp_kmax}
    Real sum2 =  (lnComp -= max).exp().dot(*p_prop_);
    // compute likelihood of each sample for each component
    p_tik_->row(i) = (*p_prop_ * lnComp.exp())/sum2;
    // compute lnLikelihood
    sum += max + std::log(sum2);
  }
  setLnLikelihood(sum);
}
/* estimate the proportions and the parameters of the components of the
 *  model given the current tik/zi mixture parameters values.
 **/
void IMixtureModelBase::mStep()
{ computeProportions();
  /* implement specific parameters estimation in concrete class. */
}

/* Compute prop using the ML estimator, default implementation. */
void IMixtureModelBase::computeProportions()
{
  for (int k=p_prop_->firstIdx(); k<= p_prop_->lastIdx(); k++)
  { p_prop_->elt(k) = p_tik_->col(k).mean();}
}

/* Compute Zi using the Map estimator, default implementation. */
void IMixtureModelBase::mapStep()
{
  for (int i = p_zi_->firstIdx(); i<= p_zi_->lastIdx(); ++i)
  {
    int k;
    p_tik_->row(i).maxElt(k);
    p_zi_->elt(i) = k;
  }
}

/* set the parameters of the  mixture model.
 *  @param p_prop pointer on the proportion of the mixture model
 *  @param p_tik pointer on the posterior probabilities
 *  @param p_zi pointer on the class labels
 * */
void IMixtureModelBase::setMixtureParameters( CArrayPoint<Real> const* p_prop
                                                  , Array2D<Real> const* p_tik
                                                  , CArrayVector<int> const* p_zi)
{
  deleteMixtureParameters();
  p_prop_ = const_cast<CArrayPoint<Real>*>(p_prop);
  p_tik_ = const_cast<Array2D<Real>*>(p_tik);
  p_zi_ = const_cast<CArrayVector<int>*>(p_zi);
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

/** delete  the mixture model parameters. */
void IMixtureModelBase::deleteMixtureParameters()
{
  if (isParametersCreated_)
  {
    if (p_prop_) delete p_prop_;  p_prop_ = 0;
    if (p_tik_) delete p_tik_;  p_tik_ = 0;
    if (p_zi_) delete p_zi_;  p_zi_ = 0;
    isParametersCreated_ = false;
  }
}

/* create the proportions */
void IMixtureModelBase::createProp()
{
  if (p_prop_)
  { p_prop_->resize(nbCluster_);}
  else
  { p_prop_ = new CArrayPoint<Real>(nbCluster_);}
  p_prop_->setValue(1./(Real)nbCluster_);
}
/* create the tik probabilities */
void IMixtureModelBase::createTik()
{
  if (p_tik_)
  { p_tik_->resize(nbSample(), nbCluster_);}
  else
  { p_tik_ = new Array2D<Real>(nbSample(), nbCluster_);}
  p_tik_->setValue(1./(Real)nbCluster_);
}
/* create the zi labels */
void IMixtureModelBase::createZi()
{
  if (p_zi_)
  { p_zi_->resize(nbSample());}
  else
  { p_zi_ = new CArrayVector<int>(nbSample());}
  p_zi_->setValue(STKBASEARRAYS);
}

} // namespace SDTK
