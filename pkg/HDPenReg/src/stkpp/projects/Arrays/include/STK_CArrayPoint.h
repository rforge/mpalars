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
 * Project:  stkpp::Arrays
 * created on: 25 nov. 2011
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_CArrayPoint.h
 *  @brief In this file we implement the final class CArrayPoint.
 **/

#ifndef STK_CARRAY2DPOINT_H
#define STK_CARRAY2DPOINT_H

#include "../../STKernel/include/STK_Constants.h"
#include "../../STKernel/include/STK_MetaTemplate.h"
#include "../../STKernel/include/STK_StaticAssert.h"

#include "../../Sdk/include/STK_Traits.h"

#include "STK_Arrays_Util.h"
#include "STK_CAllocator.h"
#include "STK_ICArray.h"

namespace STK
{
template< typename Type, int SizeCols_=UnknownSize, bool Orient_ = Arrays::by_row_>
class CArrayPoint;

template< typename Type, int SizeRows_, int SizeCols_, bool Orient_>
class CArray;
template< typename Type, int Size_, bool Orient_>
class CArraySquare;
template< typename Type, int SizeRows_, int SizeCols_, bool Orient_>
class CArrayVector;
template< typename Type, int SizeRows_, int SizeCols_, bool Orient_>
class CArrayNumber;

typedef CArrayPoint<Real, UnknownSize, Arrays::by_row_>   CPointX;
typedef CArrayPoint<Real, 2, Arrays::by_row_>             CPoint2;
typedef CArrayPoint<Real, 3, Arrays::by_row_>             CPoint3;

namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of the Traits class for CArray class.
 */
template<typename Type_, int SizeCols_, bool Orient_>
struct Traits< CArrayPoint<Type_, SizeCols_, Orient_> >
{
  private:
    class Void { };

    typedef CArrayPoint<Type_, SizeCols_, Arrays::by_col_> RowIndirect;
    typedef CArrayPoint<Type_, SizeCols_, Arrays::by_row_> RowDirect;

    typedef CArrayPoint<Type_, UnknownSize, Arrays::by_col_> SubRowIndirect;
    typedef CArrayPoint<Type_, UnknownSize, Arrays::by_row_> SubRowDirect;

    typedef CArrayVector<Type_, 1, 1, Arrays::by_row_> ColIndirect;
    typedef CArrayVector<Type_, 1, 1, Arrays::by_col_> ColDirect;

    typedef CArrayVector<Type_, 1, 1, Arrays::by_row_> SubColIndirect;
    typedef CArrayVector<Type_, 1, 1, Arrays::by_col_> SubColDirect;

    typedef CArrayPoint<Type_, UnknownSize, Arrays::by_row_> FixedRowArrayIndirect;
    typedef CArrayPoint<Type_, SizeCols_, Arrays::by_row_> FixedColArrayDirect;

    typedef CArrayPoint<Type_, UnknownSize, Arrays::by_col_> FixedRowArrayDirect;
    typedef CArrayPoint<Type_, SizeCols_, Arrays::by_col_> FixedColArrayIndirect;

  public:
    typedef CArrayNumber<Type_, 1, 1, Orient_> Number;

    typedef typename If<Orient_, RowIndirect, RowDirect >::Result  Row;
    typedef typename If<Orient_, ColDirect, ColIndirect >::Result  Col;

    typedef typename If<Orient_, SubRowIndirect, SubRowDirect >::Result  SubRow;
    typedef typename If<Orient_, SubColDirect, SubColIndirect >::Result  SubCol;

    /* Type or array (1,1) ? */
    typedef typename If<(SizeCols_ == 1), SubCol, SubRow>::Result SubVector;

    // FIXME does not seem optimal if we want only to get a subset of rows (columns)
    typedef CArrayPoint<Type_, UnknownSize, Orient_> SubArray;
//    typedef typename If< Orient_ == Arrays::by_col_
//                       , typename If<SizeRows_ != UnknownSize, FixedRowArrayDirect, FixedColArrayIndirect>::Result
//                       , typename If<SizeCols_ != UnknownSize, FixedRowArrayIndirect, FixedColArrayDirect>::Result
//                       >::Result SubArray;
    // Transposed type
    typedef CArrayVector< Type_, SizeCols_, 1, !Orient_> Transposed;
    // The CAllocator have to have the same structure than the CArray
    typedef CAllocator<Type_, Arrays::point_, 1, SizeCols_, Orient_> Allocator;

    typedef Type_ Type;

    enum
    {
      structure_ = Arrays::point_,
      orient_    = Orient_,
      sizeRows_  = 1,
      sizeCols_  = SizeCols_,
      storage_   = Arrays::dense_
    };
};

} // namespace hidden


/** @ingroup Arrays
 * @brief declaration of  the point case.
 */
template <typename Type, int SizeCols_, bool Orient_>
class CArrayPoint : public ICArray < CArrayPoint<Type, SizeCols_, Orient_> >
{
  public:
    typedef ICArray < CArrayPoint<Type, SizeCols_, Orient_> > Base;
    typedef ArrayBase < CArrayPoint<Type, SizeCols_, Orient_> > LowBase;
    enum
    {
      structure_ = Arrays::point_,
      orient_    = Orient_,
      sizeRows_  = 1,
      sizeCols_  = SizeCols_,
      storage_   = Arrays::dense_
    };

    /** Default constructor. */
    inline CArrayPoint() : Base() {}
    /** constructor with specified dimension.
     *  @param sizeCols range of the columns
     **/
    inline CArrayPoint( int const& sizeCols) : Base(1, sizeCols) {}
    /** constructor with rbeg, rend, cbeg and cend specified,
     *  initialization with a constant.
     *  @param sizeCols range of the columns
     *  @param v initial value of the container
     **/
    inline CArrayPoint( int const& sizeCols, Type const& v) : Base(1, sizeCols, v) {}
    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if T is wrapped
     **/
    inline CArrayPoint( const CArrayPoint &T, bool ref=false) : Base(T, ref) {}
    /** constructor by reference.
     *  @param allocator the allocator to wrap
     **/
    template<class OtherAllocator>
    inline CArrayPoint( OtherAllocator const& allocator) : Base(allocator) {}
    /** destructor. */
    inline ~CArrayPoint() {}
    /** operator= : set the container to a constant value.
     *  @param v the value to set
     **/
    inline CArrayPoint& operator=(Type const& v) { this->setValue(v); return *this;}
    /** operator = : overwrite the CArray with the Right hand side T.
     *  @param T the container to copy
     **/
    template<class Rhs>
    inline CArrayPoint& operator=(Rhs const& T) { return LowBase::assign(T);}
    /** operator = : overwrite the CArray with the Right hand side rhs.
     *  @param rhs the container to copy
     **/
    inline CArrayPoint& operator=(CArrayPoint const& rhs) { return LowBase::assign(rhs);}
};

} // namespace STK


#endif /* STK_CARRAY2DPOINT_H */
