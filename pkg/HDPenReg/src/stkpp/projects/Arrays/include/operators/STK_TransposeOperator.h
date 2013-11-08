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
 * Project:  stkpp::Deriveds
 * created on: 17 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_TransposeOperator.h
 *  @brief In this file we implement the TransposeOperator class.
 **/

#ifndef STK_TRANSPOSEOPERATOR_H
#define STK_TRANSPOSEOPERATOR_H

namespace STK
{

namespace hidden
{
/** @ingroup hidden
  * Trait class for transposed array.
  **/
template<int Structure_> struct TransposeTraits;

/** specialization for array2D_ */
template<> struct TransposeTraits<Arrays::array2D_>
{ enum { structure_ = Arrays::array2D_}; };
/** specialization for square_ */
template<> struct TransposeTraits<Arrays::square_>
{ enum { structure_ = Arrays::square_}; };
/** specialization for lower_triangular_ */
template<> struct TransposeTraits<Arrays::lower_triangular_>
{ enum { structure_ = Arrays::upper_triangular_}; };
/** specialization for upper_triangular_ */
template<> struct TransposeTraits<Arrays::upper_triangular_>
{ enum { structure_ = Arrays::lower_triangular_}; };
/** specialization for diagonal_ */
template<> struct TransposeTraits<Arrays::diagonal_>
{ enum { structure_ = Arrays::diagonal_}; };
/** specialization for vector_ */
template<> struct TransposeTraits<Arrays::vector_>
{ enum { structure_ = Arrays::point_}; };
/** specialization for point_ */
template<> struct TransposeTraits<Arrays::point_>
{ enum { structure_ = Arrays::vector_}; };
/** specialization for number_ */
template<> struct TransposeTraits<Arrays::number_>
{ enum { structure_ = Arrays::number_}; };

} // namespace hidden

// forward declaration
template< typename Array> class TransposeOperator;

namespace hidden
{

/** @ingroup hidden
 *  @brief Traits class for the transposed operator
 */
template<typename Derived>
struct Traits< TransposeOperator <Derived> >
{
  typedef typename Derived::Type Type;
  enum
  {
    structure_ = TransposeTraits<Derived::structure_>::structure_,
    orient_    = !Derived::orient_,
    sizeRows_  = Derived::sizeCols_,
    sizeCols_  = Derived::sizeRows_,
    storage_   = Derived::storage_
  };
};

} // end namespace hidden

// forward declaration
template<typename Derived> class TransposeOperatorBase;


/** @class TransposeOperator
  * @ingroup Arrays
  *
  * \brief Generic expression when an expression is transposed
  *
  * @tparam Array the type of the expression to which we are applying the
  * transpose operator.
  *
  * This class represents an expression where a transpose operator is applied to
  * an expression. It is the return type of the transpose operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name TransposeOperator type explicitly.
  */
template< typename Derived>
class TransposeOperator  : public TransposeOperatorBase< Derived >, public TRef<1>
{
  public:
    typedef TransposeOperatorBase< Derived > Base;
    typedef typename hidden::Traits< TransposeOperator<Derived> >::Type Type;
    enum
    {
        structure_ = hidden::Traits< TransposeOperator<Derived> >::structure_,
        orient_    = hidden::Traits< TransposeOperator<Derived> >::orient_,
        sizeRows_  = hidden::Traits< TransposeOperator<Derived> >::sizeRows_,
        sizeCols_  = hidden::Traits< TransposeOperator<Derived> >::sizeCols_,
        storage_   = hidden::Traits< TransposeOperator<Derived> >::storage_
    };
    /** Constructor */
    inline TransposeOperator( Derived const& rhs) : Base(), rhs_(rhs) {}
    /**  @return the range of the rows */
    inline Range const rows() const { return rhs_.cols();}
    /** @return the range of the Columns */
    inline Range const cols() const { return rhs_.rows();}
    /** @return the number of rows of the transposed expression */
    inline int const sizeRowsImpl() const { return rhs_.sizeCols();}
    /** @return the number of columns of the transposed expression */
    inline int const sizeColsImpl() const { return rhs_.sizeRows();}

    /** @return the right hand side expression */
    inline Derived const& rhs() const { return rhs_; }

  protected:
    Derived const& rhs_;
};

/** @ingroup Arrays
  * @brief implement the access to the elements in the (2D) general case.
  **/
template< typename Derived>
class TransposeOperatorBase : public ArrayBase< TransposeOperator< Derived> >
{
  public:
  typedef typename hidden::Traits< TransposeOperator<Derived> >::Type Type;
  enum
  {
      structure_ = hidden::Traits< TransposeOperator<Derived> >::structure_,
      orient_    = hidden::Traits< TransposeOperator<Derived> >::orient_,
      sizeRows_  = hidden::Traits< TransposeOperator<Derived> >::sizeRows_,
      sizeCols_  = hidden::Traits< TransposeOperator<Derived> >::sizeCols_,
      storage_   = hidden::Traits< TransposeOperator<Derived> >::storage_
  };
    typedef ArrayBase< TransposeOperator< Derived> > Base;
    /** constructor. */
    inline TransposeOperatorBase() : Base() {}
    /** @return the element (i,j) of the transposed expression.
     *  @param i index of the row
     *  @param j index of the column
     **/
    inline Type const elt2Impl(int i, int j) const
    { return (this->asDerived().rhs().elt(j, i));}
    /** @return the element ith element of the transposed expression
     *  @param i index of the ith element
     **/
    inline Type const elt1Impl(int i) const
    { return (this->asDerived().rhs().elt(i));}
    /** accesses to the element of the transposed expression */
    inline Type const elt0Impl() const
    { return (this->asDerived().rhs().elt());}
};

} // namespace STK

#endif /* STK_TRANSPOSEOPERATOR_H */
