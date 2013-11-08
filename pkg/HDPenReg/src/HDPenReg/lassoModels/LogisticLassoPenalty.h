/*--------------------------------------------------------------------*/
/*     Copyright (C) 2013-2013  Serge Iovleff, Quentin Grimonprez

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : quentin.grimonprez@inria.fr
*/

/*
 * Project:  MPAGenomics::
 * created on: 4 nov. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file LogisticLassoPenalty.h
 *  @brief In this file .
 **/


#ifndef LOGISTICLASSOPENALTY_H_
#define LOGISTICLASSOPENALTY_H_

#include "PenalizedModels.h"
#include "LassoPenalty.h"

namespace HD
{

  /** @ingroup lassoModels
   *  @brief The class LassoPenalty derived from the @c IPenalty class.
   *  It contains the matrix penalty associated to the lasso problem.
   */
  class LogisticLassoPenalty : public LassoPenalty
  {
    public:
      /**default constructor*/
      LogisticLassoPenalty();
      /** Constructor
       *  @param lambda penalization parameter for the l1-norm of the estimates
       *  @param n size of sample
       *  @param p size of Penalty (number of covariates)
       */
      LogisticLassoPenalty(STK::Real lambda);

      /** Copy constructor
       *  @param penalty LassoPenalty object to copy
       */
      LogisticLassoPenalty(LogisticLassoPenalty const& penalty);

      /** destructor */
      inline virtual ~LogisticLassoPenalty() {};

      /**clone*/
      LogisticLassoPenalty* clone() const;

      //getter
      /**@return A constant pointer on z*/
      inline STK::CVectorX const*  p_z() const { return &z_;}
      /**@return  z*/
      inline STK::CVectorX const&  z() const { return z_;}
      /**@return p_y_*/
      inline STK::CVectorX const*  p_y() const { return p_y_;}
      /**@return p_currentData_*/
      inline STK::CArrayXX const*  p_currentData() const { return p_currentData_;}

      //setter
      /** change the value of p_y_ */
      inline void setPy(STK::CVectorX const* p_y) {p_y_ = p_y; z_.resize(p_y->sizeRows());}
      /** change the value of p_currentData_ */
      inline void setPcurrentData(STK::CArrayXX const* p_currentData) {p_currentData_ = p_currentData;}

      //methods
      /**
       * update the lasso penalty (for fixed sigma)
       * @param beta current estimates
       */
      void update(STK::CVectorX const& beta);



    protected:
      void updateZ(STK::CVectorX const& beta);

    private:
      ///estimated response z
      STK::CVectorX z_;
      ///pointer to y
      STK::CVectorX const* p_y_;
      ///pointer to data
      STK::CArrayXX const* p_currentData_;

  };
}


#endif /* LOGISTICLASSOPENALTY_H_ */
