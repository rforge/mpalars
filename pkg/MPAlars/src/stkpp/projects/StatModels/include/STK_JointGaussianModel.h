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
 * Purpose: define the class IUnivStatModel.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_JointGaussianModel.h
 *  @brief In this file we define the class JointGaussianModel.
 **/

#ifndef STK_JOINTGAUSSIANMODEL_H
#define STK_JOINTGAUSSIANMODEL_H

#include <cmath>

#include "STK_IMultiStatModel.h"
#include "../../STatistiK/include/STK_Law_Normal.h"

namespace STK
{

/** @ingroup StatModels
 *  Structure encapsulating the parameters of a Joint Gaussian model.
 */
struct JointGaussianParameters: public IMultiParameters
{
  public:
    /** default constructor */
    JointGaussianParameters() : mu_(), sigma_() {}
    /** default constructor */
    JointGaussianParameters(Range const& range) : mu_(range, 0.), sigma_(range, 1.) {}
    /** copy constructor. @param param the parameters to copy. */
    JointGaussianParameters( JointGaussianParameters const& param)
                            : mu_(param.mu_)
                            , sigma_(param.sigma_)
    {}
    /** destructor */
    ~JointGaussianParameters() {}
    /** clone pattern. */
    JointGaussianParameters* clone() const
    { return new JointGaussianParameters(*this);}

    /** @return the means */
    inline Array2DPoint<Real> const& mu() const { return mu_;}
    /** @return the mean of the jth law */
    inline Array2DPoint<Real> const& sigma() const { return sigma_;}
    /** @return the mean of the jth law */
    inline Real const& mu(int const& j) const { return mu_[j];}
    /** @return the standard deviation of the jth law */
    inline Real const& sigma(int const& j) const { return sigma_[j];}
    /** set the mean of the jth law */
    inline void setMu(int const& j, Real const& mu) { mu_[j] = mu;}
    /** set the standard deviation of the jth law */
    inline void setSigma(int const& j, Real const& sigma) { sigma_[j] = sigma;}
    /** resize the set of parameter */
    inline void resize(Range const& size)
    { mu_.resize(size); mu_ = 0.;
      sigma_.resize(size); sigma_ = 1.;
    }
  protected:
    Array2DPoint<Real> mu_;
    Array2DPoint<Real> sigma_;
};

/** @ingroup StatModels
 * A joint Gaussian model is a statistical model of the form:
 * following form
 * \f[
 *     f(\mathbf{x}_i = x|\theta) =
 *     \prod_{j=1}^p p_{j}^{x_i^j} (1-p_{j})^{1-x_i^j},
 *      \quad x_i^j\in\{0,1\}, \quad j=1,\ldots,p, \quad i=1,\ldots,n.
 * \f]
 *
 **/
template <class Array, class WColVector = CVectorX>
class JointGaussianModel : public IMultiStatModel<Array, WColVector, JointGaussianParameters >
{
  public:
    /** Type of the data contained in the container */
    typedef typename Array::Type Type;
    /** Type of the row vector of the container */
    typedef typename Array::Row RowVector;
    /** Type of the column vector of the container */
    typedef typename Array::Col ColVector;
    /** Base class */
    typedef IMultiStatModel<Array, WColVector, JointGaussianParameters > Base;
    using Base::p_data;
    using Base::p_param;
    /** default constructor. */
    JointGaussianModel() : Base() {}
    /** Constructor with data set. */
    JointGaussianModel(Array const& data) : Base(data) {}
    /** Constructor with a ptr on the data set. */
    JointGaussianModel(Array const* p_data) : Base(p_data) {}
    /** Copy constructor. */
    JointGaussianModel(JointGaussianModel const& model) : Base(model) {}
    /** destructor */
    virtual ~JointGaussianModel(){}
    /** clone pattern. @return a clone of this. */
    JointGaussianModel* clone() const { return new JointGaussianModel(*this);}

    /** compute the number of free parameters */
    virtual int computeNbFreeParameters() const
    { return 2*p_data()->sizeCols();}
    /** compute the log Likelihood of an observation. */
    virtual Real computeLnLikelihood( RowVector const& rowData) const
    {
      Real sum =0.;
      for (Integer j= rowData.firstIdx(); j <= rowData.lastIdx(); ++j)
      {
        sum += Law::Normal::lpdf(rowData[j], p_param()->mu(j), p_param()->sigma(j));
      }
      return sum;
    }
  protected:
    /** initialize the parameters */
    virtual void initParameters()
    {
       if(!p_param())
          this->p_param_ = new JointGaussianParameters(p_data()->cols());
//       p_param()->resize(p_data()->cols());
    }
    /** compute the parameters */
    virtual void computeParameters()
    {
      for (int j=p_data()->firstIdxCols(); j<=p_data()->lastIdxCols(); ++j)
      {
        p_param()->setMu(j,Stat::mean(p_data()->col(j)));
        p_param()->setSigma(j,std::sqrt((double)Stat::varianceWithFixedMean(p_data()->col(j), p_param()->mu(j), true)));
      }
    }
    /** compute the weighted parameters */
    virtual void computeParameters( WColVector const& weights)
    {
      for (int j=p_data()->firstIdxCols(); j<=p_data()->lastIdxCols(); ++j)
      {
        p_param()->setMu(j,Stat::mean(p_data()->col(j), weights));
        p_param()->setSigma(j,std::sqrt((double)Stat::varianceWithFixedMean(p_data()->col(j), weights, p_param()->mu(j), true)));
      }
    }
};

} // namespace STK

#endif /* STK_JOINTGAUSSIANMODEL_H */
