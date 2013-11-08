/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2013 Vincent KUBICKI

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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project:  stk++
 * created on: Oct 24, 2013
 * Author:   Vincent KUBICKI
 **/

/** @file STK_Gaussian_s.h
 *  @brief In this file we implement the Gaussian_p_s and Gaussian_p_s classes
 **/

#ifndef STK_GAUSSIAN_S_H
#define STK_GAUSSIAN_S_H

#include "../STK_IMixtureModel.h"
#include "STK_DiagGaussianComponent.h"
//#include "STK_Gaussian_sImpl.h"

namespace STK
{

//forward declaration, to allow for recursive template
template<class Array>class Gaussian_pk_s;
template<class Array>class Gaussian_p_s;

namespace hidden
{
/** @ingroup hidden
 *  Traits class for the Gaussian_p_s traits policy. */
template<class _Array>
struct MixtureTraits< Gaussian_pk_s<_Array> >
{
  typedef _Array Array;
  typedef DiagGaussianComponent<_Array, Gaussian_s_Parameters> Component;
  typedef Gaussian_s_Parameters        Parameters;
};
/** @ingroup hidden
 *  Traits class for the Gaussian_p_s traits policy. */
template<class _Array>
struct MixtureTraits< Gaussian_p_s<_Array> >
{
  typedef _Array Array;
  typedef DiagGaussianComponent<_Array, Gaussian_s_Parameters> Component;
  typedef Gaussian_s_Parameters        Parameters;
};

} // namespace hidden

/** @ingroup Clustering
 *  The diagonal Gaussian mixture model @c Gaussian_p_s is
 *  the most general diagonal Gaussian model and have a density function of the
 *  form
 * \f[
 *  f(\mathbf{x}|\theta) = \sum_{k=1}^K p_k \prod_{j=1}^d
 *    \frac{1}{\sqrt{2\pi}\sigma} \exp\left\{-\frac{(x^j-\mu^j_k)^2}{2\sigma^2}\right\}.
 * \f]
 **/
template<class Array>
class Gaussian_pk_s : public IMixtureModel<Gaussian_pk_s<Array> >
{
  public:
    typedef IMixtureModel<Gaussian_pk_s<Array> > Base;
    using Base::p_data_;
    using Base::components_;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    Gaussian_pk_s( int nbCluster) : Base(nbCluster) {}
    /** constructor. Create a model with a data set.
     *  @param data the data set to process
     *  @param nbCluster number of cluster in the model
     **/
    Gaussian_pk_s( int nbCluster, Array const& data) : Base(nbCluster, data) {}
    /** constructor with a pointer on the data set
     *  @param p_data a pointer on the data set to process
     *  @param nbCluster number of cluster in the model
     **/
    Gaussian_pk_s( int nbCluster, Array const* p_data) : Base(nbCluster, p_data) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    Gaussian_pk_s( Gaussian_pk_s const& model) : Base(model) {}
    /** destructor */
    virtual ~Gaussian_pk_s() {}
    /** Initialize the component of the model.
     *  This function have to be called prior to any used of the class.
     *  In this interface, the @c initializeModel() method call the base
     *  class IMixtureModelBase::initializeModel() and for all the
     *  components initialize the shared parameter sigma_.
     **/
    virtual void initializeModel()
    {
      scale_.resize(p_data_->cols());
      IMixtureModelBase::initializeModel();
      for (int k= components_.firstIdx(); k <= components_.lastIdx(); ++k)
      { components_[k]->p_param()->p_sigma_ = &sigma_;}
    }
    /** Write the parameters*/
    virtual void writeParameters(ostream& os) const
    {
      stk_cout << _T("lnLikelihood = ") << this->lnLikelihood() << _T("\n";);
      stk_cout << _T("proportions = ") << *this->p_prop() << _T("\n";);
      for (int k= components_.firstIdx(); k <= components_.lastIdx(); ++k)
      {
        stk_cout << _T("---> Component ") << k << _T("\n";);
        stk_cout << _T("mean_ = ") << components_[k]->p_param()->mean_;
        stk_cout << _T("sigma_ = ") << sigma_ * Const::Point<Real>(this->nbVar());
      }
    }
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()*this->nbVar();}
  protected:
    Real sigma_;
};

/** @ingroup Clustering
 *  The diagonal Gaussian mixture model @c Gaussian_p_s is
 *  the most general diagonal Gaussian model and have a density function of the
 *  form
 * \f[
 *  f(\mathbf{x}|\theta) = \sum_{k=1}^K p \prod_{j=1}^d
 *    \frac{1}{\sqrt{2\pi}\sigma} \exp\left\{-\frac{(x^j-\mu^j_k)^2}{2\sigma^2}\right\}.
 * \f]
 **/
template<class Array>
class Gaussian_p_s : public IMixtureModelFixedProp<Gaussian_p_s<Array> >
{
  public:
    typedef IMixtureModelFixedProp<Gaussian_p_s<Array> > Base;
    using Base::p_data_;
    using Base::components_;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    Gaussian_p_s( int nbCluster) : Base(nbCluster) {}
    /** constructor. Create a model with a data set.
     *  @param data the data set to process
     *  @param nbCluster number of cluster in the model
     **/
    Gaussian_p_s( int nbCluster, Array const& data) : Base(nbCluster, data) {}
    /** constructor with a pointer on the data set
     *  @param p_data a pointer on the data set to process
     *  @param nbCluster number of cluster in the model
     **/
    Gaussian_p_s( int nbCluster, Array const* p_data) : Base(nbCluster, p_data) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    Gaussian_p_s( Gaussian_p_s const& model) : Base(model) {}
    /** destructor */
    virtual ~Gaussian_p_s() {}
    /** Initialize the component of the model.
     *  This function have to be called prior to any used of the class.
     *  In this interface, the @c initializeModel() method call the base
     *  class IMixtureModelBase::initializeModel() and for all the
     *  components initialize the shared parameter sigma_.
     **/
    virtual void initializeModel()
    {
      scale_.resize(p_data_->cols());
      IMixtureModelBase::initializeModel();
      for (int k= components_.firstIdx(); k <= components_.lastIdx(); ++k)
      { components_[k]->p_param()->p_sigma_ = &sigma_;}
    }

    /** Write the parameters*/
    virtual void writeParameters(ostream& os) const
    {
      stk_cout << _T("lnLikelihood = ") << this->lnLikelihood() << _T("\n";);
      stk_cout << _T("proportions = ") << *this->p_prop() << _T("\n";);
      for (int k= components_.firstIdx(); k <= components_.lastIdx(); ++k)
      {
        stk_cout << _T("---> Component ") << k << _T("\n";);
        stk_cout << _T("mean_ = ") << components_[k]->p_param()->mean_;
        stk_cout << _T("sigma_ = ") << sigma_ * Const::Point<Real>(this->nbVar());
      }
    }
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()*this->nbVar();}
  protected:
    Real sigma_;
};

} // namespace STK

#endif /* STK_GAUSSIAN_SJK_H */
