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
 * Project:  stkpp::Arrays
 * created on: 27 sept. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ArrayBaseVisitor.h
 *  @brief In this file we define the ArrayBaseVisitor classes.
 **/

#ifndef STK_ARRAYBASEVISITOR_H
#define STK_ARRAYBASEVISITOR_H

#include "../../STKernel/include/STK_MetaTemplate.h"
#include "../../STatistiK/include/STK_RandBase.h"
#include "STK_Arrays_Util.h"
#include "./visitors/STK_VisitorsImpl.h"
#include "./visitors/STK_Visitors.h"

namespace STK
{

/** Run the visitor @a visitor to the whole coefficients of the array.
  *
  * The template parameter @a Visitor is the type of the visitor and provides
  * the following interface:
  * @code
  * struct MyVisitor {
  *   // called for all  coefficients
  *   void operator() (Type const& value, Index i, Index j);
  * };
  * @endcode
  *
  * @note visitors offer automatic unrolling for small fixed size matrix.
  *
  * @sa minElt, maxElt
  */
template<typename Derived>
template<typename Visitor>
void ArrayBase<Derived>::visit(Visitor& visitor) const
{
  typedef hidden::VisitorSelector<Visitor, Derived, structure_, sizeRows_, sizeCols_> Impl;
  Impl::run(this->asDerived(), visitor);
}

/** Apply the Visitor @a visitor to the whole coefficients of the array.
  *
  * The template parameter @a Visitor is the type of the visitor and provides
  * the following interface:
  * @code
  * struct MyVisitor {
  *   // called for all  coefficients
  *   void operator() (Type& value);
  * };
  * @endcode
  *
  * @note visitors offer automatic unrolling for small fixed size matrix.
  *
  * @sa setValue
  */
template<typename Derived>
template<typename Visitor>
void ArrayBase<Derived>::visit(Visitor& visitor)
{
  typedef hidden::VisitorSelector<Visitor, Derived, structure_, sizeRows_, sizeCols_> Impl;
  Impl::apply(this->asDerived(), visitor);
}

/*
 * start public function of the ArrayBase class using visitor
 */
/* @brief Set the value to all the Allocator */
template<typename Derived>
void ArrayBase<Derived>::setValue(Type const& value)
{
  hidden::ValueVisitor<Type> visitor(value);
  visit(visitor);
}

/* @brief Set the value one to all the Array */
template<typename Derived>
void ArrayBase<Derived>::ones()
{
  hidden::ValueVisitor<Type> visitor(Type(1));
  visit(visitor);
}

/* @brief Set the value one to all the Array */
template<typename Derived>
void ArrayBase<Derived>::zeros()
{
  hidden::ValueVisitor<Type> visitor(Type(0));
  visit(visitor);
}

template<typename Derived>
void ArrayBase<Derived>::randUnif()
{
  hidden::RandUnifVisitor<Type> visitor;
  visit(visitor);
}

template<typename Derived>
void ArrayBase<Derived>::randGauss()
{
    hidden::RandGaussVisitor<Type> visitor;
    visit(visitor);
}


template<typename Derived>
typename hidden::Traits<Derived>::Type ArrayBase<Derived>::minElt( int& row, int& col) const
{
  hidden::MinEltVisitor<Type> visitor;
  visit(visitor);
  row = visitor.row_;
  col = visitor.col_;
  return visitor.res_;
}
template<typename Derived>
typename hidden::Traits<Derived>::Type ArrayBase<Derived>::maxElt( int& row, int& col) const
{
  typedef hidden::MaxEltVisitor<Type> Visitor;
  Visitor visitor;
  visit(visitor);
  row = visitor.row_;
  col = visitor.col_;
  return visitor.res_;
}
template<typename Derived>
typename hidden::Traits<Derived>::Type ArrayBase<Derived>::minElt( int& idx) const
{
  typedef hidden::MinEltVisitor<Type> Visitor;
  Visitor visitor;
  visit(visitor);
  idx = hidden::GetIdx<Visitor, hidden::Traits<Derived>::structure_ >::idx(visitor);
  return visitor.res_;
}
template<typename Derived>
typename hidden::Traits<Derived>::Type ArrayBase<Derived>::maxElt( int& idx) const
{
  typedef hidden::MaxEltVisitor<Type> Visitor;
  Visitor visitor;
  visit(visitor);
  idx = hidden::GetIdx<Visitor, hidden::Traits<Derived>::structure_ >::idx(visitor);
  return visitor.res_;
}
template<typename Derived>
typename hidden::Traits<Derived>::Type ArrayBase<Derived>::minElt() const
{
  typedef hidden::MinEltVisitor<Type> Visitor;
  Visitor visitor;
  visit(visitor);
  return visitor.res_;
}
template<typename Derived>
typename hidden::Traits<Derived>::Type ArrayBase<Derived>::maxElt() const
{
  typedef hidden::MaxEltVisitor<Type> Visitor;
  Visitor visitor;
  visit(visitor);
  return visitor.res_;
}


template<typename Derived>
typename hidden::Traits<Derived>::Type ArrayBase<Derived>::minEltSafe( int& row, int& col) const
{
  typedef hidden::MinEltSafeVisitor<Type> Visitor;
  Visitor visitor;
  visit(visitor);
  row = visitor.row_;
  col = visitor.col_;
  return visitor.res_;
}
template<typename Derived>
typename hidden::Traits<Derived>::Type ArrayBase<Derived>::maxEltSafe( int& row, int& col) const
{
  typedef hidden::MaxEltSafeVisitor<Type> Visitor;
  Visitor visitor;
  visit(visitor);
  row = visitor.row_;
  col = visitor.col_;
  return visitor.res_;
}
template<typename Derived>
typename hidden::Traits<Derived>::Type ArrayBase<Derived>::minEltSafe( int& idx) const
{
  typedef hidden::MinEltSafeVisitor<Type> Visitor;
  Visitor visitor;
  visit(visitor);
  idx = hidden::GetIdx<Visitor, hidden::Traits<Derived>::structure_ >::idx(visitor);
  return visitor.res_;
}
template<typename Derived>
typename hidden::Traits<Derived>::Type ArrayBase<Derived>::maxEltSafe( int& idx) const
{
  typedef hidden::MaxEltSafeVisitor<Type> Visitor;
  Visitor visitor;
  visit(visitor);
  idx = hidden::GetIdx<Visitor, hidden::Traits<Derived>::structure_ >::idx(visitor);
  return visitor.res_;
}
template<typename Derived>
typename hidden::Traits<Derived>::Type ArrayBase<Derived>::minEltSafe() const
{
  typedef hidden::MinEltSafeVisitor<Type> Visitor;
  Visitor visitor;
  visit(visitor);
  return visitor.res_;
}
template<typename Derived>
typename hidden::Traits<Derived>::Type ArrayBase<Derived>::maxEltSafe() const
{
  typedef hidden::MaxEltSafeVisitor<Type> Visitor;
  Visitor visitor;
  visit(visitor);
  return visitor.res_;
}



/* sum the values of all the array */
template<typename Derived>
typename hidden::Traits<Derived>::Type ArrayBase<Derived>::sum() const
{
  hidden::SumVisitor<Type> visitor;
  visit(visitor);
  return visitor.res_;
}
/* @return the norm of this*/
template<typename Derived>
typename hidden::Traits<Derived>::Type ArrayBase<Derived>::norm() const
{ return Type(std::sqrt(norm2()));}
/* @return the square norm of this*/
template<typename Derived>
typename hidden::Traits<Derived>::Type ArrayBase<Derived>::norm2() const
{ return square().sum();}
/* sum the values of all the array */
template<typename Derived>
typename hidden::Traits<Derived>::Type ArrayBase<Derived>::mean() const
{ return (this->sizeArray() >0) ? sum()/Type(this->sizeArray()) : Arithmetic<Type>::NA();}
/* sum the values of all the array */
template<typename Derived>
typename hidden::Traits<Derived>::Type ArrayBase<Derived>::variance() const
{ return (this->sizeArray() >0) ?
   ((*this-mean()).square().sum()/Type(this->sizeArray())) : Arithmetic<Type>::NA();
}

/* sum the values of all the array */
template<typename Derived>
typename hidden::Traits<Derived>::Type ArrayBase<Derived>::sumSafe() const
{
  hidden::SumSafeVisitor<Type> visitor;
  visit(visitor);
  return visitor.res_;
}
template<typename Derived>
typename hidden::Traits<Derived>::Type ArrayBase<Derived>::normSafe() const
{ return Type(std::sqrt(norm2Safe()));}
/* @return the square norm of this*/
template<typename Derived>
typename hidden::Traits<Derived>::Type ArrayBase<Derived>::norm2Safe() const
{ return square().sumSafe();}
/* sum the values of all the array */
template<typename Derived>
typename hidden::Traits<Derived>::Type ArrayBase<Derived>::meanSafe() const
{
  int size = nbValues();
  return (size>0) ? sumSafe()/Type(size) : Arithmetic<Type>::NA();
}
template<typename Derived>
typename hidden::Traits<Derived>::Type ArrayBase<Derived>::varianceSafe() const
{
  int size = nbValues();
  return (size >0) ?
    ((*this - (sumSafe()/Type(size))).square().sumSafe()/Type(size)) : Arithmetic<Type>::NA();
}


/* @return the weighted sum of all the elements of this using a Visitor*/
template<typename Derived>
template<typename Rhs>
typename hidden::Traits<Derived>::Type  ArrayBase<Derived>::wsum(ArrayBase<Rhs> const& weights) const
{ return dot(weights);}
/* @return the norm of this*/
template<typename Derived>
template<typename Rhs>
typename hidden::Traits<Derived>::Type  ArrayBase<Derived>::wnorm(ArrayBase<Rhs> const& weights) const
{ return Type(std::sqrt(wnorm2(weights)));}
/* @return the square norm of this*/
template<typename Derived>
template<typename Rhs>
typename hidden::Traits<Derived>::Type  ArrayBase<Derived>::wnorm2(ArrayBase<Rhs> const& weights) const
{ return (square().dot(weights));}
/* @return the mean of all the elements of this using a Visitor*/
template<typename Derived>
template<typename Rhs>
typename hidden::Traits<Derived>::Type  ArrayBase<Derived>::wmean(ArrayBase<Rhs> const& weights) const
{ return (this->sizeArray() >0) ? wsum(weights)/weights.sum() : Arithmetic<Type>::NA();}
/* @return the variance of all the elements of this using a Visitor*/
template<typename Derived>
template<typename Rhs>
typename hidden::Traits<Derived>::Type ArrayBase<Derived>::wvariance(ArrayBase<Rhs> const& weights) const
{
  Type size = weights.sum();
  return (size >0) ?
    ((*this-(wsum(weights)/size)).square().wsum()/size) : Arithmetic<Type>::NA();;
}

/* @return the weighted sum of all the elements of this using a Visitor*/
template<typename Derived>
template<typename Rhs>
typename hidden::Traits<Derived>::Type  ArrayBase<Derived>::wsumSafe(ArrayBase<Rhs> const& weights) const
{ return dotSafe(weights);}
/* @return the norm of this*/
template<typename Derived>
template<typename Rhs>
typename hidden::Traits<Derived>::Type  ArrayBase<Derived>::wnormSafe(ArrayBase<Rhs> const& weights) const
{ return Type(std::sqrt(wnorm2Safe(weights)));}
/* @return the square norm of this*/
template<typename Derived>
template<typename Rhs>
typename hidden::Traits<Derived>::Type  ArrayBase<Derived>::wnorm2Safe(ArrayBase<Rhs> const& weights) const
{ return (square().dotSafe(weights));}
/* @return the mean of all the elements of this using a Visitor*/
template<typename Derived>
template<typename Rhs>
typename hidden::Traits<Derived>::Type  ArrayBase<Derived>::wmeanSafe(ArrayBase<Rhs> const& weights) const
{
  Type size = weights.sum();
  return (size > 0) ? wsum(weights)/size : Arithmetic<Type>::NA();
}
/* @return the variance of all the elements of this using a Visitor*/
template<typename Derived>
template<typename Rhs>
typename hidden::Traits<Derived>::Type ArrayBase<Derived>::wvarianceSafe(ArrayBase<Rhs> const& weights) const
{
  Type size = weights.sum();
  return (size > 0) ?
    ((*this-(wsumSafe(weights)/size)).square().wsumSafe(weights)/size) : Arithmetic<Type>::NA();
}


/* count the number values in the array */
template<typename Derived>
int ArrayBase<Derived>::nbValues() const
{ return isFinite().template cast<int>().sum();}


/** @ingroup Arrays
 *  @brief Applies the visitor \a Visitor to the whole elements of the
 *  matrix or vector but don't modify the Array.
  * The template parameter \a Visitor is the type of the visitor and provides
  * the following interface:
  * @code
  * struct MyVisitor
  * {
  *   // called for all elements
  *   void operator() (const Type& value, int i, int j);
  * };
  * @endcode
  *
  * @note compared to one or two @c for loops, visitors offer automatic
  * unrolling for small fixed size matrix.
  */
template< typename Derived, typename Visitor>
class ConstArrayBaseVisitor
{
  private:
    Derived const& array_;
  public:
    typedef typename hidden::VisitorSelector< Visitor, Derived
                                        , (Arrays::Structure)hidden::Traits<Derived>::structure_
                                        , hidden::Traits<Derived>::sizeRows_
                                        , hidden::Traits<Derived>::sizeCols_> Impl;

    ConstArrayBaseVisitor( Derived const& T) : array_(T.asDerived()) { }
    ~ConstArrayBaseVisitor() { }
    inline void visit(Visitor& funct) const { Impl::run(array_, funct);}
};

/** @ingroup Arrays
 *  @brief Applies the visitor \a visitor to the whole elements of the matrix or
 *  vector.
  *
  * The template parameter \a Visitor is the type of the visitor and provides
  * the following interface:
  * @code
  * struct MyVisitor {
  *   // called for all elements
  *   void operator() (Type& value, int i, int j);
  * };
  * struct MyVisitor1D {
  *   // called for all elements
  *   void operator() (Type& value, int j);
  * };
  * @endcode
  * The value is modified by the Visitor.
  *
  * @note compared to one or two @c for loops, visitors offer automatic
  * unrolling for small fixed size matrix.
  */
template< typename Derived, typename Visitor>
class ArrayBaseVisitor
{
  private:
    Derived& array_;
  public:
    typedef typename hidden::VisitorSelector< Visitor, Derived
                                        , (Arrays::Structure)hidden::Traits<Derived>::structure_
                                        , hidden::Traits<Derived>::sizeRows_
                                        , hidden::Traits<Derived>::sizeCols_> Impl;

    ArrayBaseVisitor( Derived& T) : array_(T) { }
    ~ArrayBaseVisitor() { }
    inline void visit(Visitor& funct) { Impl::apply(array_, funct);}
};

/* @return the minimum of all elements of matrix
  * and puts in row and col its location.
  *
  * @sa  maxElt(Derived const&, int&, int&), minElt()
  */
template<typename Derived>
typename hidden::Traits<Derived>::Type minElt(Derived const& matrix, int& row, int& col)
{
  typedef typename hidden::Traits<Derived>::Type Type;
  hidden::MinEltVisitor<Type> visitor;
  ConstArrayBaseVisitor<Derived, hidden::MinEltVisitor<Type> > arrayVisitor(matrix);
  arrayVisitor.visit(visitor);
  row = visitor.row_;
  col = visitor.col_;
  return (visitor.res_);
}

/* @returns the maximum of all elements of matrix
  * and puts in row and col its location.
  *
  * @sa ArrayBase::minElt(int), ArrayBase::maxElt(int,int), ArrayBase::visitor(), ArrayBase::minElt()
  */
template<typename Derived>
typename hidden::Traits<Derived>::Type maxElt(Derived const& matrix, int& row, int& col)
{
  typedef typename hidden::Traits<Derived>::Type Type;
  hidden::MaxEltVisitor<Type> visitor;
  ConstArrayBaseVisitor<Derived, hidden::MaxEltVisitor<Type> > arrayVisitor(matrix);
  arrayVisitor.visit(visitor);
  row = visitor.row_;
  col = visitor.col_;
  return visitor.res_;
}

/* @return the minimum of all elements of *this
  * and puts in index its location.
  *
  * @sa minElt(int,int), ArrayBase::maxElt(int,int), ArrayBase::visitor(), ArrayBase::minElt()
  */
template<typename Derived>
typename hidden::Traits<Derived>::Type minElt(Derived const& tab, int& idx)
{
  typedef typename hidden::Traits<Derived>::Type Type;
  typedef hidden::MinEltVisitor<Type> Visitor;
  Visitor visitor;
  ConstArrayBaseVisitor<Derived, hidden::MinEltVisitor<Type> > arrayVisitor(tab);
  arrayVisitor.visit(visitor);
  idx = hidden::GetIdx<Visitor, hidden::Traits<Derived>::structure_ >::idx(visitor);
  return visitor.res_;
}

/* @return the maximum of all elements of tab
  * and puts in index its location.
  *
  * @sa ArrayBase::maxElt(int&,int&), ArrayBase::minElt(int&,int&)
  * , ArrayBase::visitor(), ArrayBase::maxElt()
  */
template<typename Derived>
typename hidden::Traits<Derived>::Type maxElt(Derived const& tab, int& idx)
{
  typedef typename hidden::Traits<Derived>::Type Type;
  /* assert nn row or nb column is 1*/
  typedef hidden::MaxEltVisitor<Type> Visitor;
  Visitor visitor;
  ConstArrayBaseVisitor<Derived, hidden::MaxEltVisitor<Type> > arrayVisitor(tab);
  arrayVisitor.visit(visitor);
  idx = hidden::GetIdx<Visitor, hidden::Traits<Derived>::structure_ >::idx(visitor);
  return visitor.res_;
}

} // namespace STK

#endif /* STK_ARRAYBASEVISITOR_H */
