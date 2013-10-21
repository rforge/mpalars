/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2013  Serge Iovleff

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

    Contact : Serge.Iovleff@stkpp.org
*/

/*
 * Project: stkpp::Clustering
 * created on: 5 sept. 2013
 * Author:  iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_Gamma_ajk_bj.h
 *  @brief In this file we define the Gamma_pk_ajk_bj and Gamma_p_ajk_bj mixture classes.
 **/

#ifndef STK_GAMMA_AJK_BJ_H
#define STK_GAMMA_AJK_BJ_H

#include "../STK_IMixtureModel.h"

#include "STK_GammaComponent.h"
#include "STK_Gamma_ajk_bjImpl.h"

namespace STK
{
template<class Array>class Gamma_pk_ajk_bj;
template<class Array>class Gamma_p_ajk_bj;

namespace hidden
{
/** @ingroup hidden
 *  Traits class for the Gamma_pk_ajk_bj traits policy. */
template<class _Array>
struct MixtureTraits< Gamma_pk_ajk_bj<_Array> >
{
  typedef _Array Array;
  typedef GammaComponent<_Array, Gamma_ajk_bj_Parameters> Component;
  typedef Gamma_ajk_bj_Parameters        Parameters;
};

/** @ingroup hidden
 *  Traits class for the Gamma_p_ajk_bjk traits policy. */
template<class _Array>
struct MixtureTraits< Gamma_p_ajk_bj<_Array> >
{
  typedef _Array Array;
  typedef GammaComponent<_Array, Gamma_ajk_bj_Parameters> Component;
  typedef Gamma_ajk_bj_Parameters        Parameters;
};

} // namespace hidden

/** @ingroup Clustering
 *  Gamma_pk_ajk_bj is a mixture model of the following form
 * \f[
 *     f(\mathbf{x}_i|\theta) = \sum_{k=1}^K p_k
 *     \prod_{j=1}^p\left(\frac{x_i^j}{b_{j}}\right)^{a_{jk}-1}
 *                   \frac{e^{-x_i^j/b_{j}}} {b_{j} \, \Gamma(a_{jk})},
 *      \quad x_i^j>0, \quad i=1,\ldots,n.
 * \f]
 **/
template<class Array>
class Gamma_pk_ajk_bj : public IMixtureModel<Gamma_pk_ajk_bj<Array> >
{
  public:
    typedef IMixtureModel<Gamma_pk_ajk_bj<Array> > Base;
    using Base::p_data_;
    using Base::components_;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    Gamma_pk_ajk_bj( int nbCluster) : Base(nbCluster), scale_() {}
    /** constructor. Create a model with a data set.
     *  @param data the data set to process
     *  @param nbCluster number of cluster in the model
     **/
    Gamma_pk_ajk_bj( int nbCluster, Array const& data)
                     : Base(nbCluster, data), scale_(data.cols())
    {}
    /** constructor with a pointer on the data set
     *  @param p_data a pointer on the data set to process
     *  @param nbCluster number of cluster in the model
     **/
    Gamma_pk_ajk_bj( int nbCluster, Array const* p_data) : Base(nbCluster, p_data)
    {}
    /** copy constructor
     *  @param model The model to copy
     **/
    Gamma_pk_ajk_bj( Gamma_pk_ajk_bj const& model) : Base(model), scale_(model.scale_)
    {}
    /** destructor */
    virtual ~Gamma_pk_ajk_bj() {}
    /** Initialize the component of the model.
     *  This function have to be called prior to any used of the class.
     *  In this interface, the @c initializeModel() method call the base
     *  class IMixtureModelBase::initializeModel() and for all the
     *  components create the parameters.
     **/
    virtual void initializeModel()
    {
      scale_.resize(p_data_->cols());
      IMixtureModelBase::initializeModel();
      for (int k= components_.firstIdx(); k <= components_.lastIdx(); ++k)
      { components_[k]->p_param()->p_scale_ = &scale_;}
    }
    /** Write the parameters*/
    virtual void writeParameters(ostream& os) const
    {
      stk_cout << "lnLikelihood = " << this->lnLikelihood() << "\n";
      stk_cout << "proportions = " << *this->p_prop() << "\n";
      for (int k= components_.firstIdx(); k <= components_.lastIdx(); ++k)
      {
        stk_cout << "---> Component " << k << "\n";
        stk_cout << "shape_ = "
                 << components_[k]->p_param()->shape_;
        stk_cout << "scale_ = "
                 << scale_;
      }
    }
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()*this->nbVar()+ this->nbVar() + this->nbCluster()-1;}

  protected:
    /** Array of the common scale */
    Array2DPoint<Real> scale_;
};

/** @ingroup Clustering
 *  Gamma_p_ajk_bj is a mixture model of the following form
 * \f[
 *     f(\mathbf{x}_i|\theta) = \sum_{k=1}^K \frac{1}{K}
 *     \prod_{j=1}^p\left(\frac{x_i^j}{b_{j}}\right)^{a_{jk}-1}
 *                   \frac{e^{-x_i^j/b_{j}}} {b_{j} \, \Gamma(a_{jk})},
 *      \quad x_i^j>0, \quad j=1,\ldots,p, \quad i=1,\ldots,n.
 * \f]
 **/
template<class Array>
class Gamma_p_ajk_bj : public IMixtureModelFixedProp< Gamma_p_ajk_bj<Array> >
{
  public:
    typedef IMixtureModelFixedProp< Gamma_p_ajk_bj<Array> > Base;
    using Base::p_data_;
    using Base::components_;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    Gamma_p_ajk_bj( int nbCluster) : Base(nbCluster)
    {}
    /** constructor. Create a model with a data set.
     *  @param data the data set to process
     *  @param nbCluster number of cluster in the model
     **/
    Gamma_p_ajk_bj( int nbCluster, Array const& data) : Base(nbCluster, data), scale_(data.cols())
    {}
    /** constructor with a pointer on the data set
     *  @param p_data a pointer on the data set to process
     *  @param nbCluster number of cluster in the model
     **/
    Gamma_p_ajk_bj( int nbCluster, Array const* p_data) : Base(nbCluster, p_data)
    {}
    /** copy constructor
     *  @param model The model to copy
     **/
    Gamma_p_ajk_bj( Gamma_p_ajk_bj const& model) : Base(model), scale_(model.scale_)
    {}
    /** destructor */
    virtual ~Gamma_p_ajk_bj() {}
    /** Initialize the component of the model.
     *  This function have to be called prior to any used of the class.
     *  In this interface, the @c initializeModel() method call the base
     *  class IMixtureModelBase::initializeModel() and for all the
     *  components create the parameters.
     **/
    virtual void initializeModel()
    {
      scale_.resize(p_data_->cols());
      IMixtureModelBase::initializeModel();
      for (int k= components_.firstIdx(); k <= components_.lastIdx(); ++k)
      { components_[k]->p_param()->p_scale_ = &scale_;}
    }
    /** Write the parameters*/
    virtual void writeParameters(ostream& os) const
    {
      stk_cout << "lnLikelihood = " << this->lnLikelihood() << "\n";
      stk_cout << "proportions = " << *this->p_prop() << "\n";
      for (int k= components_.firstIdx(); k <= components_.lastIdx(); ++k)
      {
        stk_cout << "---> Component " << k << "\n";
        stk_cout << "shape_ = "
                 << components_[k]->p_param()->shape_;
        stk_cout << "scale_ = "
                 << scale_;
      }
    }
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()*this->nbVar()+ this->nbVar();}

  protected:
    /** Array of the common scale */
    Array2DPoint<Real> scale_;
};

}  // namespace STK

#endif /* STK_GAMMA_AJK_BJ_H */
