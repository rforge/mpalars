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
 * created on: 13 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ArrayBase.h
 *  @brief In this file we define the base class for Arrays and Expressions
 **/

#ifndef STK_ARRAYBASE_H
#define STK_ARRAYBASE_H

/// utility macro allowing to construct binary operators
#define MAKE_BINARY_OPERATOR(OPERATOR, FUNCTOR) \
  template<typename Rhs> \
  inline BinaryOperator< FUNCTOR<Type, typename hidden::Traits<Rhs>::Type>, Derived, Rhs> \
  (OPERATOR)(const ArrayBase<Rhs> &other) const \
  { return BinaryOperator<FUNCTOR<Type, typename hidden::Traits<Rhs>::Type>, Derived, Rhs>(this->asDerived(), other.asDerived()) ;}


// forward declarations
namespace STK
{
template<class Derived> class ArrayBase;
template<class Derived, class Rhs> struct  ProductReturnType;
template<class Derived> class  ArrayInitializer;
} // namespace STK

#include "../../Sdk/include/STK_Traits.h"

#include "../../STKernel/include/STK_MetaTemplate.h"
#include "../../STKernel/include/STK_Functors.h"
#include "../../STKernel/include/STK_Macros.h"
#include "../../STKernel/include/STK_StaticAssert.h"
#include "../../STKernel/include/STK_Exceptions.h"
#include "../../STKernel/include/STK_Range.h"

#include "STK_Arrays_Util.h"
#include "STK_ITContainer2D.h"
#include "STK_IContainerRef.h"

#include "./products/STK_ProductOperators.h"

#include "./operators/STK_TransposeOperator.h"
#include "./operators/STK_UnaryOperators.h"
#include "./operators/STK_BinaryOperators.h"
#include "./operators/STK_DotOperators.h"

namespace STK
{

/** @ingroup Arrays
 *  @brief base class for templated evaluation expressions and visitors.
 *
 * This class is the base that is inherited by all objects (matrix, vector,
 * point and expression). The common API for these objects is contained in
 * this class.
 *
 *  @tparam Derived is the derived type, e.g., a matrix, vector, point type or
 *  an expression.
 **/

template<class Derived>
class ArrayBase : public ITContainer<Derived>
{
  public:
    typedef ITContainer<Derived> Base;
    typedef typename hidden::Traits<Derived>::Type Type;

    enum
    {
      structure_ = hidden::Traits<Derived>::structure_,
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };

  protected:
    /** Default constructor */
    inline ArrayBase() : Base() {}
    /** destructor */
    inline ~ArrayBase() {}

  public:
    /** Visit the container using a constant visitor. @sa ConstArrayBaseVisitor
     *  @param visitor the visitor to run
     **/
    template<typename Visitor>
    void visit(Visitor& visitor) const;
    /** Visit the container using a visitor. @sa ArrayBaseVisitor
     *  @param visitor the visitor to apply
     **/
    template<typename Visitor>
    void visit(Visitor& visitor);
    /** set a value to this container. @sa visit(), ones(), zeros()
     *  @param value the value to set
     **/
    void setValue(Type const& value);
    /** set random values to this using a uniform law. @sa visit(), ones(), zeros() */
    void randUnif();
    /** set random values to this using a standard gaussian law. @sa visit(), ones(), zeros() */
    void randGauss();
    /** set one to this using a Visitor. @sa visit(), setValue(), zeros() */
    void ones();
    /** set zero to this using a Visitor. @sa visit(), ones(), setValue()*/
    void zeros();

    /** @return the number of available values in the array (not count NA values).*/
    int nbValues() const;

    /** @return the minimum of all elements of this using a Visitor
      * and puts in (row, col) its location.
      * @sa maxElt(int,int), visitor(), minElt()
      */
    Type minElt( int& row, int& col) const;
    /** @return the maximum of all elements of this using a Visitor
      * and puts in (row, col) its location.
      * @sa minElt(), visitor()
      */
    Type maxElt( int& row, int& col) const;
    /** @return the minimum of all elements of this using a Visitor
      * and puts in  pos its location.
      * @note Have to be used for col-vector, row-vector or diagonal matrix only.
      * @sa maxElt(), visitor()
      */
    Type minElt( int& pos) const;
    /** @return the maximum of all elements of this using a Visitor
      * and puts in pos its location.
      * @note Have to be used for col-vector, row-vector or diagonal matrix only.
      * @sa minElt(), visitor()
      */
    Type maxElt( int& pos) const;
    /** @return the minimum of all elements of this using a Visitor.
      * @sa maxElt(row, col), visitor()
      */
    Type minElt() const;
    /** @return the maximum of all elements of this using a Visitor.
      * @sa minElt(row, col), visitor()
      */
    Type maxElt() const;

    /** @return the minimum of all elements of this using a Visitor
      * and puts in (row, col) its location.
      * @sa maxElt(int,int), visitor(), minElt()
      */
    Type minEltSafe( int& row, int& col) const;
    /** @return the maximum of all elements of this using a Visitor
      * and puts in (row, col) its location.
      * @sa minElt(), visitor()
      */
    Type maxEltSafe( int& row, int& col) const;
    /** @return the minimum of all elements of this using a Visitor
      * and puts in  pos its location.
      * @note Have to be used for col-vector, row-vector or diagonal matrix only.
      * @sa maxElt(), visitor()
      */
    Type minEltSafe( int& pos) const;
    /** @return the maximum of all elements of this using a Visitor
      * and puts in pos its location.
      * @note Have to be used for col-vector, row-vector or diagonal matrix only.
      * @sa minElt(), visitor()
      */
    Type maxEltSafe( int& pos) const;
    /** @return the minimum of all elements of this using a Visitor
      * @sa maxElt(row, pos), visitor()
      */
    Type minEltSafe() const;
    /** @return the maximum of all elements of this using a Visitor
      * @sa minElt(row, pos), visitor()
      */
    Type maxEltSafe() const;

    /** @return the sum of all the elements of this using a Visitor*/
    inline Type sum() const;
    /** @return the norm of this*/
    inline Type norm() const;
    /** @return the square norm of this*/
    inline Type norm2() const;
    /** @return the mean of all the elements of this using a Visitor*/
    inline Type mean() const;
    /** @return the variance of all the elements of this using a Visitor*/
    inline Type variance() const;

    /** @return the sum of all the elements of this using a Visitor*/
    inline Type sumSafe() const;
    /** @return the norm of this*/
    inline Type normSafe() const;
    /** @return the square norm of this*/
    inline Type norm2Safe() const;
    /** @return the mean of all the elements of this using a safe Visitor*/
    inline Type meanSafe() const;
    /** @return the variance of all the elements of this using a safe Visitor*/
    inline Type varianceSafe() const;

    /** @return the weighted sum of all the elements of this using a Visitor
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    inline Type wsum(ArrayBase<Rhs> const& weights) const;
    /** @return the weighted norm of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    inline Type wnorm(ArrayBase<Rhs> const& weights) const;
    /** @return the weighted square norm of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    inline Type wnorm2(ArrayBase<Rhs> const& weights) const;
    /** @return the weighted mean of all the elements of this using a Visitor
     * @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    inline Type wmean(ArrayBase<Rhs> const& weights) const;
    /** @return the weighted variance of all the elements of this using a Visitor
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    inline Type wvariance(ArrayBase<Rhs> const& weights) const;

    /** @return the weighted sum of all the elements of this using a Visitor
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    inline Type wsumSafe(ArrayBase<Rhs> const& weights) const;
    /** @return the weighted norm of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    inline Type wnormSafe(ArrayBase<Rhs> const& weights) const;
    /** @return the weighted square norm of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    inline Type wnorm2Safe(ArrayBase<Rhs> const& weights) const;
    /** @return the weighted mean of all the elements of this using a Visitor
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    inline Type wmeanSafe(ArrayBase<Rhs> const& weights) const;
    /** @return the weighted variance of all the elements of this using a Visitor
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    inline Type wvarianceSafe(ArrayBase<Rhs> const& weights) const;

    /** @return a copy of @c rhs inside @c this object.
     *  If the ranges of @c this and @c rhs are not exactly the same, the assign
     *  method will call the resize method on this.
     *
     *  @note If @c this is a reference, it cannot be resized and thus an
     *  exception will be thrown.
     **/
    template<class Rhs> Derived& assign(ArrayBase<Rhs> const& rhs);

    /** @return the matrix or vector obtained by setting this constant*/
    inline Derived& operator=( Type const& rhs)
    { this->setValue(rhs); return this->asDerived();}
    /** @return the matrix or vector obtained by evaluating this expression */
    inline Derived& operator=( Derived const& rhs)
    { return this->assign(rhs);}
    /** @return the matrix or vector obtained by evaluating this expression */
    template<typename Rhs>
    inline Derived& operator=( ArrayBase<Rhs> const& rhs)
    { return this->assign(rhs.asDerived());}
    /** Adding a Rhs to this. */
    template<typename Rhs>
    Derived& operator+=( ArrayBase<Rhs> const& other);
    /** substract a Rhs to this. */
    template<typename Rhs>
    Derived& operator-=( ArrayBase<Rhs> const& other);
    /** Adding a constant to this. */
    Derived& operator+=( Type const other);
    /** Substract a constant to this. */
    Derived& operator-=( Type const other);
    /** product of this by a constant. */
    Derived& operator*=( Type const other);
    /** dividing this by a constant. */
    Derived& operator/=( Type const other);

    /** @return an expression from the difference of this and  other. */
    MAKE_BINARY_OPERATOR(operator-,DifferenceOp)
    /** @return an expression from the addition operator. */
    MAKE_BINARY_OPERATOR(operator+,SumOp)
    /** @return an expression from the quotient operator. */
    MAKE_BINARY_OPERATOR(operator/,DivOp)
    /** @return an expression from min of this and  other */
    MAKE_BINARY_OPERATOR(min,MinOp)
    /** @return an expression from the max of this and  other */
    MAKE_BINARY_OPERATOR(max,MaxOp)
    /** @return an expression from the \< operator of *this and  other*/
    MAKE_BINARY_OPERATOR(operator<,LessOp)
    /** @return an expression from the \<= operator of *this and  other */
    MAKE_BINARY_OPERATOR(operator<=,LeqOp)
    /** @return an expression from the \> operator of *this and  other */
    MAKE_BINARY_OPERATOR(operator>,GreaterOp)
    /** @return an expression from the \>= operator of *this and  other */
    MAKE_BINARY_OPERATOR(operator>=,GeqOp)
    /** @return an expression from the == operator of *this and  other */
    MAKE_BINARY_OPERATOR(operator==,EqualOp)
    /** @return an expression from the != operator of *this and  other */
    MAKE_BINARY_OPERATOR(operator!=,NotEqualOp)

    /** @return an expression of this with each elements incremented by
     *  the constant number */
    inline UnaryOperator<AddOp<Type>, Derived>
    operator+(Type const number) const
    { return UnaryOperator<AddOp<Type>, Derived>(this->asDerived(), AddOp<Type>(number));}
    // handle the case number + array
    friend inline UnaryOperator<AddOp<Type>, Derived>
    operator+(Type const number, ArrayBase<Derived> const& other)
    { return other.asDerived() + number;}
    /** @return an expression of this with each elements decremented by the
     * constant  number */
    inline UnaryOperator<AddOp<Type>, Derived>
    operator-(Type const number) const
    { return UnaryOperator<AddOp<Type>, Derived>(this->asDerived(), AddOp<Type>(-number));}
    // handle the case number - array
    friend inline UnaryOperator<AddOp<Type>, UnaryOperator<OppositeOp<Type>, Derived> >
    operator-(Type const number, const ArrayBase<Derived>& other)
    { return (-other.asDerived()) + number;}
    /** @return an expression of the opposite of this */
    inline UnaryOperator<OppositeOp<Type>, Derived>
    operator-() const { return UnaryOperator<OppositeOp<Type>, Derived>(this->asDerived()); }
    /** @return an expression of this scaled by the number factor number */
    inline UnaryOperator< MultipleOp<Type>, Derived>
    operator*(Type const number) const
    { return UnaryOperator<MultipleOp<Type>, Derived>(this->asDerived(), MultipleOp<Type>(number));}
    // handle the case number * array
    inline friend UnaryOperator< MultipleOp<Type>, Derived>
    operator*(Type const number, ArrayBase<Derived> const& matrix)
    { return matrix.asDerived()*number; }
    /** @return an expression of this divided by the number value number */
    inline UnaryOperator<QuotientOp<Type>, Derived>
    operator/(Type const number) const
    { return UnaryOperator<QuotientOp<Type>, Derived>(this->asDerived(), QuotientOp<Type>(number));}
    /** @return which values of this is a NA value */
    inline UnaryOperator<IsNaOp<Type>, Derived> isNA() const
    { return UnaryOperator<IsNaOp<Type>, Derived>(this->asDerived());}
    /** @return which values of this is a NA value */
    inline UnaryOperator<NotOp<Type>, Derived> Not() const
    { return UnaryOperator<NotOp<Type>, Derived>(this->asDerived());}
    /** @return which values of this are finite value */
    inline UnaryOperator<IsFiniteOp<Type>, Derived> isFinite() const
    { return UnaryOperator<IsFiniteOp<Type>, Derived>(this->asDerived()); }
    /** @return which values of this are finite value */
    inline UnaryOperator<IsInfiniteOp<Type>, Derived> isInfinite() const
    { return UnaryOperator<IsInfiniteOp<Type>, Derived>(this->asDerived()); }
    /** @return an expression of the absolute value of this */
    inline UnaryOperator<AbsOp<Type>, Derived> abs() const
    { return UnaryOperator<AbsOp<Type>, Derived>(this->asDerived()); }
    /** @return an expression of the exponential of this. */
    inline UnaryOperator<ExpOp<Type>, Derived> exp() const
    { return UnaryOperator<ExpOp<Type>, Derived>(this->asDerived());}
    /** @return an expression of the logarithm of this. */
    inline UnaryOperator<LogOp<Type>, Derived> log() const
    { return UnaryOperator<LogOp<Type>, Derived>(this->asDerived());}
    /** @return an expression of the square root of this. */
    inline UnaryOperator<SqrtOp<Type>, Derived> sqrt() const
    { return UnaryOperator<SqrtOp<Type>, Derived>(this->asDerived());}
    /** @return an expression of the cosine of this. */
    inline UnaryOperator<CosOp<Type>, Derived> cos() const
    { return UnaryOperator<CosOp<Type>, Derived>(this->asDerived());}
    /** @return an expression of the sine of this.  */
    inline UnaryOperator<SinOp<Type>, Derived> sin() const
    { return UnaryOperator<SinOp<Type>, Derived>(this->asDerived());}
    /** @return an expression of the arc cosine of this. */
    inline UnaryOperator<AcosOp<Type>, Derived> acos() const
    { return UnaryOperator<AcosOp<Type>, Derived>(this->asDerived());}
    /** @return an expression of the arc sine of this. */
    inline UnaryOperator<AsinOp<Type>, Derived> asin() const
    { return UnaryOperator<AsinOp<Type>, Derived>(this->asDerived());}
    /** @return an expression of the tan of this. */
    inline UnaryOperator<TanOp<Type>, Derived> tan() const
    { return UnaryOperator<TanOp<Type>, Derived> (this->asDerived());}
    /** @return an expression of the power of this to the given exponent. */
    inline UnaryOperator<PowOp<Type>, Derived> pow( Type const exponent) const
    { return UnaryOperator<PowOp<Type>, Derived> (this->asDerived(), PowOp<Type>(exponent));}
    /** @return an expression of the inverse of this.  */
    inline UnaryOperator<InverseOp<Type>, Derived> inverse() const
    { return UnaryOperator<InverseOp<Type>, Derived>(this->asDerived());}
    /** @return an expression of the square of this. */
    inline UnaryOperator<SquareOp<Type>, Derived> square() const
    { return UnaryOperator<SquareOp<Type>, Derived>(this->asDerived());}
    /** @return an expression of the cube of this. */
    inline UnaryOperator<CubeOp<Type>, Derived> cube() const
    { return UnaryOperator<CubeOp<Type>, Derived>(this->asDerived());}
    /** @return an expression of *this == number. */
    inline UnaryOperator<EqualOp<Type, Type>, Derived> operator==(Type const number) const
    { return UnaryOperator<EqualOp<Type, Type>, Derived>(this->asDerived(), EqualOp<Type, Type>(number));}
    /** @return an expression of *this != number. */
    inline UnaryOperator<NotEqualOp<Type, Type>, Derived> operator!=(Type const number) const
    { return UnaryOperator<NotEqualOp<Type, Type>, Derived> (this->asDerived(), NotEqualOp<Type, Type>(number));}
    /** @return an expression of *this < number. */
    inline UnaryOperator<LessOp<Type, Type>, Derived> operator<(Type const number) const
    { return UnaryOperator<LessOp<Type, Type>, Derived>(this->asDerived(), LessOp<Type, Type>(number));}
    /** @return an expression of *this <= number. */
    inline UnaryOperator<LeqOp<Type, Type>, Derived> operator<=(Type const number) const
    { return UnaryOperator<LeqOp<Type, Type>, Derived>(this->asDerived(), LeqOp<Type, Type>(number));}
    /** @return an expression of *this > number. */
    inline UnaryOperator<GreaterOp<Type, Type>, Derived> operator>(Type const number) const
    { return UnaryOperator<GreaterOp<Type, Type>, Derived>(this->asDerived(), GreaterOp<Type, Type>(number));}
    /** @return the expression of *this >= number. */
    inline UnaryOperator<GeqOp<Type, Type>, Derived> operator>=(Type const number) const
    { return UnaryOperator<GeqOp<Type, Type>, Derived> (this->asDerived(), GeqOp<Type, Type>(number));}
    /** @return an expression of *this with the  Type type casted to  OtherType. */
    template<typename CastedType>
    inline UnaryOperator<CastOp<Type, CastedType>, Derived> cast() const
    { return UnaryOperator<CastOp<Type, CastedType>, Derived>(this->asDerived());}
    /** @return an expression of funct0 to this. */
    template< template<typename> class OtherOperator>
    inline UnaryOperator<OtherOperator<Type>, Derived> funct0() const
    { return UnaryOperator<OtherOperator<Type>, Derived>(this->asDerived());}
    /** @return an expression of funct1 to this. */
    template< template<typename, typename> class OtherOperator>
    inline UnaryOperator<OtherOperator<Type, Type>, Derived> funct1(Type const number) const
    { return UnaryOperator<OtherOperator<Type, Type>, Derived>(this->asDerived(), OtherOperator<Type, Type>(number));}
    /** @return the expression of transpose of this. */
    inline TransposeOperator<Derived> transpose() const
    { return TransposeOperator<Derived> (this->asDerived());}

    /** @returns the dot product of this with other.
      * @sa squaredNorm(), norm(), DotProduct
      */
    template<class Rhs>
    typename hidden::Promote<Type, typename Rhs::Type>::result_type const
    dot(ArrayBase<Rhs> const& other) const;
    /** @returns the safe dot product of this with other.
      * @sa squaredNorm(), norm(), DotProduct
      */
    template<class Rhs>
    typename hidden::Promote<Type, typename Rhs::Type>::result_type const
    dotSafe(ArrayBase<Rhs> const& other) const;

    /** @return the matricial product of this with other.*/
    template<typename Rhs>
    typename ProductReturnType<Derived, Rhs>::ReturnType const
    operator*( ArrayBase<Rhs> const& other) const;

    /** Convenient operator to set the coefficients of a matrix.
      *
      * The coefficients must be provided in the row/column order and exactly
      * match the size of the matrix. Otherwise an exception is throwed.
      */
    ArrayInitializer<Derived> operator<<(Type const& s);

    /** \sa operator<<(Type const&) */
    template<typename Rhs>
    ArrayInitializer<Derived> operator<<(ArrayBase<Rhs> const& other);
};

/** @ingroup Arrays
 *  @brief Interface class for homogeneous 2D containers which can be
 *  part of an expression and handle their own storage (like IArray2D derived
 *  classes).
 *  @sa IArray2DBase, Const::Identity
 **/
template < class Derived, int SizeRow_ = hidden::Traits<Derived>::sizeRows_
                        , int SizeCol_ = hidden::Traits<Derived>::sizeCols_>
class ITArrayBase : protected TContainer2D<SizeRow_, SizeCol_>
                  , public ArrayBase<Derived>
{
  protected:
    /** Type of the Base container */
    typedef TContainer2D<SizeRow_, SizeCol_ > Base2D;
    /** Type of the Base container */
    typedef ArrayBase<Derived> Base;
    /** Default constructor. Default values are cols=(1:0) and rows=(1:0). */
    inline ITArrayBase() : Base2D(), Base() {}
    /** constructor with specified Range.
     *  @param I the vertical range
     *  @param J the horizontal range
     **/
    inline ITArrayBase( Range const& I, Range const& J) : Base2D(I, J), Base() {}
    /** Copy constructor.
     *  @param T the container to copy
     **/
    inline ITArrayBase( ITArrayBase const& T) : Base2D(T), Base() {}
    /** destructor. */
    inline ~ITArrayBase() {}

  public:
    /**@return the Horizontal range */
    inline Range cols() const { return Base2D::cols();}
    /** @return the index of the first column */
    inline int const& firstIdxCols() const { return Base2D::firstIdxCols();}
    /**  @return the index of the last column */
    inline int lastIdxCols() const { return Base2D::lastIdxCols();}
    /** @return the Horizontal size (the number of column) */
    inline int sizeColsImpl() const { return Base2D ::sizeColsImpl();}
    /** @return the Vertical range */
    inline Range rows() const { return Base2D::rows();}
    /** @return the index of the first row */
    inline int firstIdxRows() const { return Base2D::firstIdxRows();}
    /** @return the index of the last row */
    inline int lastIdxRows() const { return Base2D::lastIdxRows();}
    /** @return the Vertical size (the number of rows) */
    inline int sizeRowsImpl() const { return Base2D::sizeRowsImpl();}
    /**  @return @c true if the container is empty, @c false otherwise */
    inline bool empty() const { return Base2D::empty();}
};

#undef MAKE_BINARY_OPERATOR

} // namespace STK

#endif /* STK_ARRAYBASE_H_ */
