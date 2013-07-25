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

/** @file STK_JointGammaModel.h
 *  @brief In this file we define the class JointGammaModel.
 **/

#ifndef STK_JOINTGAMMAMODEL_H
#define STK_JOINTGAMMAMODEL_H

#include <cmath>

#include "STK_IMultiStatModel.h"
#include "../../STatistiK/include/STK_Law_Gamma.h"
#include "../../Analysis/include/STK_Algo_FindZero.h"
#include "../../Analysis/include/STK_Funct_raw.h"

namespace STK
{

/** @ingroup StatModels
 *  Structure encapsulating the parameters of a Joint Gamma model.
 */
struct JointGammaParameters: public IMultiParameters<JointGammaParameters>
{
  public:
    /** default constructor */
    JointGammaParameters() : a_(), b_() {}
    /** default constructor */
    JointGammaParameters(Range const& range) : a_(range, 1.), b_(range, 1.) {}
    /** copy constructor. @param param the parameters to copy. */
    JointGammaParameters( JointGammaParameters const& param)
                        : a_(param.a_)
                        , b_(param.b_)
    {}
    /** destructor */
    ~JointGammaParameters() {}

    /** @return the means */
    inline Array2DPoint<Real> const& shape() const { return a_;}
    /** @return the mean of the jth law */
    inline Array2DPoint<Real> const& scale() const { return b_;}
    /** @return the mean of the jth law */
    inline Real const& shape(int const& j) const { return a_[j];}
    /** @return the standard deviation of the jth law */
    inline Real const& scale(int const& j) const { return b_[j];}
    /** set the mean of the jth law */
    inline void setShape(int const& j, Real const& shape) { a_[j] = shape;}
    /** set the standard deviation of the jth law */
    inline void setScale(int const& j, Real const& scale) { b_[j] = scale;}
    /** resize the set of parameter */
    inline void resizeImpl(Range const& size)
    { a_.resize(size); a_ = 0.;
      b_.resize(size); b_ = 1.;
    }
    /** print the parameters a_ and b_.
     *  @param os the output stream for the parameters
     **/
    inline void printImpl(ostream &os)
     { os << a_ << b_ << _T("\n");}

  protected:
    Array2DPoint<Real> a_;
    Array2DPoint<Real> b_;
};

/** @ingroup StatModels
 * A joint Gamma model is a statistical model of the
 * following form
 * \f[
 *     f(\mathbf{x}_i|\theta) =
 *     \prod_{j=1}^p\left(\frac{x_i^j}{b_j}\right)^{a_j-1}
 *                   \frac{e^{-x_i^j/b_j}}{b_j \, \Gamma(a_j)},
 *      \quad x_i^j\in\{0,1\}, \quad j=1,\ldots,p, \quad i=1,\ldots,n.
 * \f]
 **/

template <class Array, class WColVector = CVectorX>
class JointGammaModel : public IMultiStatModel<Array, WColVector, JointGammaParameters >
{

  public:
    /** Type of the data contained in the container */
    typedef typename Array::Type Type;
    /** Type of the row vector of the container */
    typedef typename Array::Row RowVector;
    /** Type of the column vector of the container */
    typedef typename Array::Col ColVector;
    /** Base class */
    typedef IMultiStatModel<Array, WColVector, JointGammaParameters > Base;
    using Base::p_data;
    using Base::p_param;
    /** default constructor. */
    JointGammaModel() : Base() {}
    /** Constructor with data set. */
    JointGammaModel(Array const& data) : Base(data){}
    /** Constructor with a ptr on the data set. */
    JointGammaModel(Array const* p_data) : Base(p_data){}
    /** Copy constructor. */
    JointGammaModel(JointGammaModel const& model) : Base(model) {}
    /** destructor */
    virtual ~JointGammaModel(){}
    /** clone pattern. @return a clone of this. */
    JointGammaModel* clone() const { return new JointGammaModel(*this);}

    /** @return the vector of the mean of the observations */
    RowVector const& mean() const {return mean_;}
    /** vector of the mean log of the observations */
    RowVector const& meanLog() const {return meanLog_;}
    /** vector of the variance of the observations */
    RowVector const& variance() const {return variance_;}
    /** compute the number of free parameters */
    virtual int computeNbFreeParameters() const
    { return 2*p_data()->sizeCols();}
    /** compute the log Likelihood of an observation. */
    virtual Real computeLnLikelihood( RowVector const& rowData) const
    {
      Real sum =0.;
      for (Integer j= rowData.firstIdx(); j <= rowData.lastIdx(); ++j)
      { sum += Law::Gamma::lpdf(rowData[j], p_param()->shape(j), p_param()->scale(j));}
      return sum;
    }
  protected:
    class dloglikelihood : public IFunction<dloglikelihood >
    {
      public:
        dloglikelihood( Real const& mean, Real const& meanLog)
                      : delta_(meanLog - Real(std::log(mean))) {}
        /** @return the value of the function at a
         * @param a a positive real value
         **/
        inline Real fImpl(Real a) const
        { return (delta_ + std::log(a) - Funct::psi_raw(a));}
        /** @return the minimal value of the function at x */
        inline Real xminImpl() const { return 0;}
      private:
        Real delta_;
    };
    /** compute the parameters */
    virtual void computeParameters()
    {
      computeMeans();
      for (int j=p_data()->firstIdxCols(); j<=p_data()->lastIdxCols(); ++j)
      {
        Real start1 = (mean_[j]*mean_[j]) / variance_[j];
        Real start2 = 0.9*start1 +  0.05/(mean_[j] - meanLog_[j]);
        dloglikelihood funct(mean_[j], meanLog_[j]);
        Real shape =  Algo::findZero(funct, start1, start2);
        // replace with moment estimator if needed
        if (!Arithmetic<Real>::isFinite(shape)) { shape =  mean_[j]*mean_[j]/variance_[j];}
        p_param()->setShape(j, shape);
        p_param()->setScale(j, mean_[j]/shape);
      }
    }
    /** compute the weighted parameters */
    virtual void computeParameters( WColVector const& weights)
    {
      computeMeans(weights);
      for (int j=p_data()->firstIdxCols(); j<=p_data()->lastIdxCols(); ++j)
      {
        Real start1 = (mean_[j]*mean_[j]) / variance_[j];
        Real start2 = 0.9*start1 +  0.05/(mean_[j] - meanLog_[j]);
        dloglikelihood funct(mean_[j], meanLog_[j]);
        Real shape =  Algo::findZero(funct, start1, start2);
        // replace with moment estimator if needed
        if (!Arithmetic<Real>::isFinite(shape)) { shape =  mean_[j]*mean_[j]/variance_[j];}
        p_param()->setShape(j, shape);
        p_param()->setScale(j, mean_[j]/shape);
      }
    }

  private:
    /** vector of the mean of the observations */
    RowVector mean_;
    /** vector of the mean log of the observations */
    RowVector meanLog_;
    /** vector of the variance of the observations */
    RowVector variance_;
    /** compute the mean and the mean log of the ith observations */
    void computeMeans()
    {
      mean_.resize(p_data()->cols());
      meanLog_.resize(p_data()->cols());
      variance_.resize(p_data()->cols());
      for (int j=p_data()->firstIdxCols(); j<=p_data()->lastIdxCols(); ++j)
      {
        mean_[j] =  p_data()->col(j).meanSafe();
        meanLog_[j] = p_data()->col(j).log().meanSafe();
        variance_[j] = p_data()->col(j).varianceSafe();
      }
    }
    /** compute the mean and the mean log of the ith observations */
    void computeMeans(WColVector const& weights)
    {
      mean_.resize(p_data()->cols());
      meanLog_.resize(p_data()->cols());
      variance_.resize(p_data()->cols());
      for (int j=p_data()->firstIdxCols(); j<=p_data()->lastIdxCols(); ++j)
      {
        mean_[j] =  p_data()->col(j).wmeanSafe(weights);
        meanLog_[j] = p_data()->col(j).log().wmeanSafe(weights);
        variance_[j] = p_data()->col(j).wvarianceSafe(weights);
      }
    }

};

} // namespace STK

#endif /* STK_JOINTGAMMAMODEL_H */
