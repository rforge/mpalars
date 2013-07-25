/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2007  Serge Iovleff

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
 * Project:  stkpp::STatistiK::StatDesc
 * Purpose:  Compute elementary 1D statistics for all variables.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Stat_Univariate.h
 *  @brief This file contain the specialization of the class Univariate
 *  for the Real Type.
 **/

#ifndef STK_STAT_FUNCTORS_H
#define STK_STAT_FUNCTORS_H

#include "../../STKernel/include/STK_MetaTemplate.h"
#include "../../Arrays/include/STK_ITContainer.h"
#include "../../DManager/include/STK_HeapSort.h"

#include "STK_Stat_Univariate.h"

namespace STK
{

namespace Stat
{

template<typename Array, bool isVector > struct MinOp;
template<typename Array, bool isVector > struct MaxOp;
template<typename Array, bool isVector > struct MeanOp;
template<typename Array, bool isVector > struct VarianceOp;
template<typename Array, bool isVector > struct VarianceWithFixedMeanOp;

} // namespace Stat

namespace hidden
{

template<typename Array, template<class, bool> class StatOp>
struct StatOpSelector
{
  enum
  {
    isVector_   =  (  Array::structure_ == int(Arrays::vector_)
                   || Array::structure_ == int(Arrays::point_)
                   )
  };
  typedef StatOp<Array, (bool)isVector_> TypeOp;
  typedef typename TypeOp::result_type result_type;
};

/** @ingroup StatDesc
 *  Compute the minimal value of the variable V discarding all missing values.
 *  @param V the variable
 *  @return the minimal value of the variable V
 **/
template<class TContainer1D>
Real minImpl( ArrayBase<TContainer1D> const&  V)
{
  Real min  = Arithmetic<Real>::infinity();
  for (int i=V.firstIdx(); i<=V.lastIdx(); i++)
  { if (!Arithmetic<Real>::isNA(V[i])) min = std::min(min, V[i]);}
  return min;
}

/** @ingroup StatDesc
 *  Compute the maximal value of the variable V discarding all missing values.
 *  @param V the variable
 *  @return the maximal value of the variable V
 **/
template<class TContainer1D>
Real maxImpl( ArrayBase<TContainer1D> const&  V)
{
  Real max  = - Arithmetic<Real>::infinity();
  for (int i=V.firstIdx(); i<=V.lastIdx(); i++)
  { if (!Arithmetic<Real>::isNA(V[i])) max = std::max(max, V[i]);}
  return max;
}

/** @ingroup StatDesc
 *  Compute the mean of the variable V discarding all missing values.
 *  \f[ \hat{\mu} = \frac{1}{n} \sum_{i=1}^n V(i) \f]
 *  @param V the variable
 *  @return the mean or NA if there is no available value
 **/
template<class TContainer1D>
Real meanImpl( ArrayBase<TContainer1D> const&  V)
{
  // no samples
  if (V.empty()) { return Arithmetic<Real>::NA();}
  // get dimensions
  int nobs = V.size();
  // sum the samples
  Real sum  = 0.0;
  for (int i=V.firstIdx(); i<=V.lastIdx(); i++)
  { (!Arithmetic<Real>::isNA(V[i])) ? sum += V[i] : nobs--;}
  // compute the mean
  return nobs ? sum /= Real(nobs) : Arithmetic<Real>::NA();
}

/** @ingroup StatDesc
 *  Compute the variance of the variable V with fixed mean.
 *  \f[ \hat{\sigma^2} = \frac{1}{n} \sum_{i=1}^n (V(i) - \mu)^2. \f]
 *  using a compensated algorithm.
 *
 *  @note
 *  Chan, Tony F.; Golub, Gene H.; LeVeque, Randall J. (1983).
 *  Algorithms for Computing the Sample Variance: Analysis and Recommendations.
 *  The American Statistician 37, 242-247.
 *
 *  @param V variable
 *  @param mu the fixed mean
 *  @param unbiased @c true if we want an unbiased estimator of the variance,
 *  @c false otherwise
 **/
template<class TContainer1D>
Real varianceWithFixedMeanImpl( ArrayBase<TContainer1D> const& V, Real const& mu, bool unbiased)
{
  // no samples
  if (V.empty()) { return Arithmetic<Real>::NA();}

  int nobs = V.size();
  // sum
  Real sum  = 0.0, var  = 0.0, dev;
  for (int i=V.firstIdx(); i<=V.lastIdx(); i++)
  {
    if (!Arithmetic<Real>::isNA(V[i]))
    {
      sum += (dev = V[i] - mu); // deviation from the mean
      var += (dev*dev);         // squared value
    }
    else nobs--;
  }
  // compute the variance
  if (unbiased)
  {
    return (nobs > 1) ? (var - (sum*sum)/Real(nobs))/Real(nobs -1)
                      : Arithmetic<Real>::NA();
  }
  // variance
  return (nobs > 0) ? (var - (sum*sum)/(Real)nobs)/(Real)(nobs)
                    : Arithmetic<Real>::NA();
}

/** @ingroup StatDesc
 *  Compute the variance of the variable V.
 *  \f[ \hat{\sigma}^2 = \frac{1}{n} \sum_{i=1}^n (V(i)-\hat{\mu})^2. \f]
 *  @param V variable
 *  @param unbiased @c true if we want an unbiased estimator of the variance,
 *  @c false otherwise
 **/
template<class TContainer1D>
Real varianceImpl( ArrayBase<TContainer1D> const& V, bool unbiased)
{
  // no samples
  if (V.empty()) { return Arithmetic<Real>::NA();}

  int nobs = V.size();
  // Compute the mean
  Real mu = meanImpl<TContainer1D>(V);
  // sum
  Real sum  = 0.0, var  = 0.0, dev;
  for (int i=V.firstIdx(); i<=V.lastIdx(); i++)
  {
    if (!Arithmetic<Real>::isNA(V[i]))
    {
      sum += (dev = V[i] - mu); // deviation from the mean
      var += (dev*dev);         // squared value
    }
    else nobs--;
  }
  // compute the variance
  if (unbiased)
  {
    return (nobs > 1) ? (var - (sum*sum)/Real(nobs))/Real(nobs -1)
                      : Arithmetic<Real>::NA();
  }
  // variance
  return (nobs > 0) ? (var - (sum*sum)/(Real)nobs)/(Real)(nobs)
                    : Arithmetic<Real>::NA();
}


/** @ingroup StatDesc
 *  Compute the (weighted) mean of the variable V
 *  \f[ \hat{\mu} = \frac{1}{\sum_{i=1}^n W(i)} \sum_{i=1}^n W(i) V(i). \f]
 *  If the range of the weights does not match the range
 *  of the variable the method return the usual mean.
 *  @param V variable
 *  @param W weights
 **/
template<class TContainer1D, class WColVector>
Real wmeanImpl( ArrayBase<TContainer1D> const& V, WColVector const& W)
{
  // no samples
  if (V.empty()||W.empty()) { return Arithmetic<Real>::NA();}
  // use common range
  Range common = Range::inf(V.range(), W.range());

  // sum the weighted samples
  Real sum  = 0.0, nweight = 0.0;
  for (int i=common.firstIdx(); i<=common.lastIdx(); i++)
  { if ( (!Arithmetic<Real>::isNA(V[i])) && (!Arithmetic<Real>::isNA(W[i])))
    {
      Real weight  = std::abs(W[i]);
      nweight  += weight;
      sum      += weight * V[i];
    }
  }
  // compute the weighted mean. If all weights are 0, we get 0
  return (nweight) ? sum /= nweight : 0.;
}

/** Compute the weighted variance of the variable V with fixed mean.
 *  \f[ \hat{\sigma^2} = \frac{1}{\sum_{i=1}^n W(i)}
 *                  \sum_{i=1}^n W(i) (V(i) - \mu)^2
 *  \f]
 *  @param V variable
 *  @param W weights
 *  @param mu the mean
 *  @param unbiased @c true if we want an unbiased estimator of the variance,
 *  @c false otherwise
 **/
template<class TContainer1D, class WColVector>
Real wvarianceWithFixedMeanImpl( ArrayBase<TContainer1D> const& V, WColVector const& W
                               , Real const& mu, bool unbiased
                               )
{
  // no samples
  if (V.empty()) { return Arithmetic<Real>::NA();}
  // use common range
  Range common = Range::inf(V.range(), W.range());
  // sum the weighted samples
  Real dev, sum = 0.0, var = 0.0, nweight = 0.0, nweight2 = 0.0;
  for (int i=common.firstIdx(); i<=common.lastIdx(); i++)
  { if ( !Arithmetic<Real>::isNA(V[i]) && !Arithmetic<Real>::isNA(W[i]) )
    {
      Real weight = std::abs(W[i]);
      nweight    += weight;
      nweight2   += weight * weight;
      sum        += weight*(dev = V[i]-mu); // deviation from the mean
      var        += weight*(dev*dev);       // squared value
    }
  }
  // compute the variance
  if (unbiased)
  {
    return (nweight*nweight - nweight2 > 0.) ? (var - sum*sum/nweight)/(nweight - nweight2/nweight)
                                             : 0.;

  }
  return (nweight) ? (var - sum*sum)/(nweight) : 0.;
}

/** @ingroup StatDesc
 *  Compute the weighted variance of the variable V.
 *  \f[ \hat{\sigma}^2
 *  = \frac{\sum_{i=1}^n W(i)}{\left( \sum_{i=1}^n W(i))\right)^2-\sum_{i=1}^n W(i)^2}
 *    \sum_{i=1}^n W(i) (V(i)-\hat{\mu})^2.
 *  \f]
 * If there is no weights, this definition reduces to the usual
 * definition of the variance with factor 1/(n-1). If the range of the weights
 * is not equal to the range of the varaible, the usual varaince is computed.
 *  @param V variable
 *  @param W weights
 *  @param unbiased @c true if we want an unbiased estimator of the variance,
 *  @c false otherwise
 **/
template<class TContainer1D, class WColVector>
Real wvarianceImpl( ArrayBase<TContainer1D> const& V, WColVector const& W
                  , bool unbiased
                  )
{
  // no samples
  if (V.empty()) { return Arithmetic<Real>::NA();}
  // Compute the mean if necessary
  Real mu = wmeanImpl(V, W);
  // use common range
  Range common = Range::inf(V.range(), W.range());
  // sum the weighted samples
  Real dev, sum = 0.0, var = 0.0, nweight = 0.0, nweight2 = 0.0;
  for (int i=common.firstIdx(); i<=common.lastIdx(); i++)
  { if ( !Arithmetic<Real>::isNA(V[i]) && !Arithmetic<Real>::isNA(W[i]) )
    {
      Real weight = std::abs(W[i]);
      nweight    += weight;
      nweight2   += weight * weight;
      sum        += weight*(dev = V[i]-mu); // deviation from the mean
      var        += weight*(dev*dev);       // squared value
    }
  }
  // compute the variance
  if (unbiased)
  {
    return (nweight*nweight - nweight2 > 0.) ? (var - sum*sum/nweight)/(nweight - nweight2/nweight)
                                             : 0.;

  }
  return (nweight) ? (var - sum*sum)/(nweight) : 0.;
}

} // namespace hidden

namespace Stat
{
/** @ingroup StatDesc
 *  @brief compute the minimal value by column.
 **/
template<typename Array>
struct MinOp<Array, false>
{
  enum { NbParam_ = 1 };
  typedef Array const& param1_type ;
  typedef typename hidden::Promote< typename Array::Type, Real >::result_type Type;
  typedef Array2DPoint<Type> result_type;

  inline MinOp( ArrayBase<Array> const& lhs) : lhs_(lhs.asDerived()) {}
  inline result_type operator()()
  {
    result_type res(lhs_.cols());
    for (int j= lhs_.firstIdxCols(); j <= lhs_.lastIdxCols(); ++j)
    { res[j] = hidden::minImpl(lhs_.col(j));}
    return res;
  }
  template<typename Weights>
  inline result_type operator()(ArrayBase<Weights> const& weights)
  {
    result_type res_(lhs_.cols());
    for (int j= lhs_.firstIdxCols(); j <= lhs_.lastIdxCols(); ++j)
    { res_[j] = hidden::minImpl(lhs_.col(j));}
    return res_;
  }
  protected:
    param1_type lhs_;
};

/** @ingroup StatDesc
 *  @brief compute the min of a row or column oriented vector.
 **/
template<typename Array>
struct MinOp<Array, true>
{
  enum { NbParam_ = 1 };
  typedef Array const& param1_type ;
  typedef typename hidden::Promote< typename Array::Type, Real >::result_type Type;
  typedef Type result_type;

  inline MinOp( ArrayBase<Array> const& lhs) : lhs_(lhs.asDerived()) {}
  inline result_type operator()() { return hidden::minImpl(lhs_);}
  template<typename Weights>
  inline result_type operator()(ArrayBase<Weights> const& weights)
  { return hidden::minImpl(lhs_);}

  protected:
    param1_type lhs_;
};

/** @ingroup StatDesc
 *  @brief compute the max by column.
 **/
template<typename Array>
struct MaxOp<Array, false>
{
  enum { NbParam_ = 1 };
  typedef Array const& param1_type ;
  typedef typename hidden::Promote< typename Array::Type, Real >::result_type Type;
  typedef Array2DPoint<Type> result_type;

  inline MaxOp( ArrayBase<Array> const& lhs) : lhs_(lhs.asDerived()) {}
  inline result_type operator()()
  {
    result_type res(lhs_.cols());
    for (int j= lhs_.firstIdxCols(); j <= lhs_.lastIdxCols(); ++j)
    { res[j] = hidden::maxImpl(lhs_.col(j));}
    return res;
  }
  template<typename Weights>
  inline result_type operator()(ArrayBase<Weights> const& weights)
  {
    result_type res_(lhs_.cols());
    for (int j= lhs_.firstIdxCols(); j <= lhs_.lastIdxCols(); ++j)
    { res_[j] = hidden::maxImpl(lhs_);}
    return res_;
  }
  protected:
    param1_type lhs_;
};

/** @ingroup StatDesc
 *  @brief compute the maximal value of a row or column oriented vector.
 **/
template<typename Array>
struct MaxOp<Array, true>
{
  enum { NbParam_ = 1 };
  typedef Array const& param1_type ;
  typedef typename hidden::Promote< typename Array::Type, Real >::result_type Type;
  typedef Type result_type;

  inline MaxOp( ArrayBase<Array> const& lhs) : lhs_(lhs.asDerived()) {}
  inline result_type operator()() { return hidden::maxImpl(lhs_);}
  template<typename Weights>
  inline result_type operator()(ArrayBase<Weights> const& weights)
  { return hidden::maxImpl(lhs_);}

  protected:
    param1_type lhs_;
};

/** @ingroup StatDesc
 *  @brief compute the mean by column.
 **/
template<typename Array>
struct MeanOp<Array, false>
{
  enum { NbParam_ = 1 };
  typedef Array const& param1_type ;
  typedef typename hidden::Promote< typename Array::Type, Real >::result_type Type;
  typedef Array2DPoint<Type> result_type;

  inline MeanOp( ArrayBase<Array> const& lhs)
               : lhs_(lhs.asDerived())
  {}
  inline result_type operator()()
  {
    result_type res(lhs_.cols());
    for (int j= lhs_.firstIdxCols(); j <= lhs_.lastIdxCols(); ++j)
    { res[j] = hidden::meanImpl(lhs_.col(j));}
    return res;
  }
  template<typename Weights>
  inline result_type operator()(ArrayBase<Weights> const& weights)
  {
    result_type res(lhs_.cols());
    for (int j= lhs_.firstIdxCols(); j <= lhs_.lastIdxCols(); ++j)
    { res[j] = hidden::wmeanImpl(lhs_.col(j), weights);}
    return res;
  }
  protected:
    param1_type lhs_;
};

/** @ingroup StatDesc
 *  @brief compute the mean by column.
 **/
template<typename Array>
struct MeanOp<Array, true>
{
  enum { NbParam_ = 1 };
  typedef Array const& param1_type ;
  typedef typename hidden::Promote< typename Array::Type, Real >::result_type Type;
  typedef Type result_type;

  inline MeanOp( ArrayBase<Array> const& lhs) : lhs_(lhs.asDerived()) {}
  inline result_type operator()() { return hidden::meanImpl(lhs_);}
  template<typename Weights>
  inline result_type operator()(ArrayBase<Weights> const& weights)
  { return hidden::wmeanImpl(lhs_, weights);}

  protected:
    param1_type lhs_;
};

/** @ingroup StatDesc
 *  @brief compute the mean by column.
 **/
template<typename Array>
struct VarianceOp<Array, false>
{
  enum { NbParam_ = 1 };
  typedef Array const& param1_type ;
  typedef typename hidden::Promote< typename Array::Type, Real >::result_type Type;
  typedef Array2DPoint<Type> result_type;

  inline VarianceOp( ArrayBase<Array> const& lhs) : lhs_(lhs.asDerived()) {}
  inline result_type operator()(bool unbiased)
  {
    result_type res(lhs_.cols());
    for (int j= lhs_.firstIdxCols(); j <= lhs_.lastIdxCols(); ++j)
    { res[j] = hidden::varianceImpl(lhs_.col(j), unbiased);}
    return res;
  }
  template<typename Weights>
  inline result_type operator()(ArrayBase<Weights> const& weights, bool unbiased)
  {
    result_type res(lhs_.cols());
    for (int j= lhs_.firstIdxCols(); j <= lhs_.lastIdxCols(); ++j)
    { res[j] = hidden::wvarianceImpl(lhs_.col(j), weights, unbiased);}
    return res;
  }
  protected:
    param1_type lhs_;
};

/** @ingroup StatDesc
 *  @brief compute the mean by column.
 **/
template<typename Array>
struct VarianceOp<Array, true>
{
  enum { NbParam_ = 1 };
  typedef Array const& param1_type ;
  typedef typename hidden::Promote< typename Array::Type, Real >::result_type Type;
  typedef Type result_type;

  inline VarianceOp( ArrayBase<Array> const& lhs) : lhs_(lhs.asDerived()) {}
  inline result_type operator()(bool unbiased) { return hidden::varianceImpl(lhs_, unbiased);}
  template<typename Weights>
  inline result_type operator()(ArrayBase<Weights> const& weights, bool unbiased)
  { return hidden::wvarianceImpl(lhs_, weights, unbiased);}

  protected:
    param1_type lhs_;
};

/** @ingroup StatDesc
 *  @brief compute the mean by column.
 **/
template<typename Array>
struct VarianceWithFixedMeanOp<Array, false>
{
  enum { NbParam_ = 1 };
  typedef Array const& param1_type ;
  typedef typename hidden::Promote< typename Array::Type, Real >::result_type Type;
  typedef typename MeanOp<Array, false>::result_type const& mean_type;
  typedef Array2DPoint<Type> result_type;

  inline VarianceWithFixedMeanOp( ArrayBase<Array> const& lhs, mean_type mean)
               : lhs_(lhs.asDerived())
               , mean_(mean)
  {}
  inline result_type operator()(bool unbiased)
  {
    result_type res(lhs_.cols());
    for (int j= lhs_.firstIdxCols(); j <= lhs_.lastIdxCols(); ++j)
    { res[j] = hidden::varianceWithFixedMeanImpl(lhs_.col(j), mean_[j], unbiased);}
    return res;
  }
  template<typename Weights>
  inline result_type operator()(ArrayBase<Weights> const& weights, bool unbiased)
  {
    result_type res(lhs_.cols());
    for (int j= lhs_.firstIdxCols(); j <= lhs_.lastIdxCols(); ++j)
    { res[j] = hidden::wvarianceWithFixedMeanImpl(lhs_.col(j), weights, mean_[j], unbiased);}
    return res;
  }
  protected:
    param1_type lhs_;
    mean_type mean_;
};

/** @ingroup StatDesc
 *  @brief compute the mean by column.
 **/
template<typename Array>
struct VarianceWithFixedMeanOp<Array, true>
{
  enum { NbParam_ = 1 };
  typedef Array const& param1_type ;
  typedef typename hidden::Promote< typename Array::Type, Real >::result_type Type;
  typedef Type result_type;

  inline VarianceWithFixedMeanOp( ArrayBase<Array> const& lhs, Type mean)
                                : lhs_(lhs.asDerived())
                                , mean_(mean)
  {}
  inline result_type operator()(bool unbiased)
  { return hidden::varianceWithFixedMeanImpl(lhs_, mean_, unbiased);}
  template<typename Weights>
  inline result_type operator()(ArrayBase<Weights> const& weights, bool unbiased)
  { return hidden::wvarianceWithFixedMeanImpl(lhs_, weights, mean_, unbiased);}

  protected:
    param1_type lhs_;
    Type mean_;
};

/** Compute the minimal(s) value(s) of A. If A is a row-vector or a
 *  column-vector then the function will return the usual min of the vector.
 *  If A is a two-dimensional array, the function will return (by value) an
 *  Array2DPoint with the minimal values of each columns.
 *  @sa mean, max, variance, varianceWithFixedMean, hidden::minImpl.
 *  @param A the array
 *  @return the minimal value(s) of A or NA if there is no available
 *  value.
 **/
template< class Array>
typename hidden::StatOpSelector<Array, MinOp>::result_type min(Array const& A)
{ return typename hidden::StatOpSelector<Array, MinOp>::TypeOp(A)();}

/** Compute the weighted minimal(s) value(s) of A. This give the same result
 *  than min(A).
 *  @sa mean, max, variance, varianceWithFixedMean, hidden::minImpl.
 *  @param A the array
 *  @param W the weights
 *  @return the minimal value(s) of A or NA if there is no available
 *  value.
 **/
template< class Array, class WColVector>
typename hidden::StatOpSelector<Array, MinOp>::result_type min(Array const& A, WColVector W)
{ return typename hidden::StatOpSelector<Array, MinOp>::TypeOp(A)(W);}

/** Compute the maximal(s) value(s) of A. If A is a row-vector or a
 *  column-vector then the function will return the usual max of the vector.
 *  If A is a two-dimensional array, the function will return (by value) an
 *  Array2DPoint with the maximal values of each columns.
 *  @sa mean, min, variance, varianceWithFixedMean, hidden::maxImpl.
 *  @param A the array
 *  @return the maximal value(s) of A or NA if there is no available
 *  value.
 **/
template< class Array>
typename hidden::StatOpSelector<Array, MaxOp>::result_type max(Array const& A)
{ return typename hidden::StatOpSelector<Array, MaxOp>::TypeOp(A)();}

/** Compute the weighted maximal(s) value(s) of A. This give the same result
 *  than max(A).
 *  @sa mean, min, variance, varianceWithFixedMean, hidden::maxImpl.
 *  @param A the array
 *  @param W the weights
 *  @return the maximal value(s) of A or NA if there is no available
 *  value.
 **/
template< class Array, class WColVector>
typename hidden::StatOpSelector<Array, MaxOp>::result_type max(Array const& A, WColVector W)
{ return typename hidden::StatOpSelector<Array, MaxOp>::TypeOp(A)(W);}

/** Compute the mean(s) value(s) of A. If A is a row-vector or a
 *  column-vector then the function will return the usual mean of the vector.
 *  If A is a two-dimensional array, the function will return (by value) an
 *  Array2DPoint with the mean values of each columns.
 *  @sa max, min, variance, varianceWithFixedMean, hidden::meanImpl.
 *  @param A the Array
 *  @return the mean(s) or NA if there is no available value
 **/
template< class Array>
typename hidden::StatOpSelector<Array, MeanOp>::result_type mean(Array const& A)
{ return typename hidden::StatOpSelector<Array, MeanOp>::TypeOp(A)();}

/** Compute the weighted mean(s) value(s) of A. This give the same result
 *  than mean(A).
 *  @sa max, min, variance, varianceWithFixedMean, hidden::meanImpl.
 *  @param A the Array
 *  @param W the weights
 *  @return the weighted mean(s) or NA if there is no available value
 **/
template< class Array, class WColVector>
typename hidden::StatOpSelector<Array, MeanOp>::result_type mean(Array const& A, WColVector W)
{ return typename hidden::StatOpSelector<Array, MeanOp>::TypeOp(A)(W);}

/** Compute the variance(s) value(s) of A. If A is a row-vector or a
 *  column-vector then the function will return the usual variance of the vector.
 *  If A is a two-dimensional array, the function will return (by value) an
 *  Array2DPoint with the variance values of each columns.
 *  @sa meax, min, mean, varianceWithFixedMean, hidden::varianceImpl.
 *  @param A the Array
 *  @param unbiased the unbiased variance(s) if @c true or the Maximum-likelihood
 *  variance otherwise (the default)
 *  @return the variance(s) or NA if there is no available value
 **/
template< class Array>
typename hidden::StatOpSelector<Array, VarianceOp>::result_type variance(Array const& A, bool unbiased = false)
{ return typename hidden::StatOpSelector<Array, VarianceOp>::TypeOp(A)(unbiased);}

/** Compute the weighted variance(s) value(s) of A. This give the same result
 *  than variance(A).
 *  @sa max, min, mean, varianceWithFixedMean, hidden::varianceImpl.
 *  @param A the Array
 *  @param W the weights
 *  @param unbiased the unbiased variance(s) if @c true or the Maximum-likelihood
 *  variance otherwise (the default)
 *  @return the variance(s) or NA if there is no available value
 **/
template< class Array, class WColVector>
typename hidden::StatOpSelector<Array, VarianceOp>::result_type variance(Array const& A, WColVector W, bool unbiased = false)
{ return typename hidden::StatOpSelector<Array, VarianceOp>::TypeOp(A)(W, unbiased);}

/** Compute the VarianceWithFixedMean(s) value(s) of A. If A is a row-vector or a
 *  column-vector then the function will return the usual variance of the vector.
 *  If A is a two-dimensional array, the function will return (by value) an
 *  Array2DPoint with the variance values of each columns.
 *  @sa meax, min, mean, varianceWithFixedMean, hidden::varianceImpl.
 *  @param A the Array
 *  @param mean The mean (s) to use
 *  @param unbiased the unbiased variance(s) if @c true or the Maximum-likelihood
 *  variance otherwise (the default)
 *  @return the variance(s) or NA if there is no available value
 **/
template< class Array>
typename hidden::StatOpSelector<Array, VarianceWithFixedMeanOp>::result_type
  varianceWithFixedMean(Array const& A, typename hidden::StatOpSelector<Array, MeanOp >::result_type mean, bool unbiased = false)
{ return typename hidden::StatOpSelector<Array, VarianceWithFixedMeanOp>::TypeOp(A, mean)(unbiased);}

/** Compute the weighted VarianceWithFixedMean(s) value(s) of A. This give the same result
 *  than mean(A).
 *  @sa max, min, mean, varianceWithFixedMean, hidden::wvarianceWithFixedMeanImpl.
 *  @param A the Array
 *  @param mean The mean (s) to use
 *  @param W the weights
 *  @param unbiased the unbiased variance(s) if @c true or the Maximum-likelihood
 *  variance otherwise (the default)
 *  @return the variance(s) or NA if there is no available value
 **/
template< class Array, class WColVector>
typename hidden::StatOpSelector<Array, VarianceWithFixedMeanOp>::result_type
  varianceWithFixedMean(Array const& A, WColVector W, typename hidden::StatOpSelector<Array, MeanOp >::result_type mean, bool unbiased = false)
{ return typename hidden::StatOpSelector<Array, VarianceWithFixedMeanOp>::TypeOp(A, mean)(W, unbiased);}

}  // namespace Stat

}  // namespace STK

#endif /*STK_STAT_FUNCTORS_H*/
