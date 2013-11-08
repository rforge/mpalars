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

/** @file LogisticLassoPenalty.cpp
 *  @brief In this file .
 **/


#include "LogisticLassoPenalty.h"
#include "../../stkpp/projects/STatistiK/include/STK_Law_Normal.h"

namespace HD
{
  /*default cosntructor*/
  LogisticLassoPenalty::LogisticLassoPenalty()
                            : LassoPenalty()
                            , z_()
                            , p_y_(0)
                            , p_currentData_(0)
  {
  }

  /* Constructor
  *  @param lambda penalization parameter for the l1-norm of the estimates
  *  @param n size of sample
  *  @param p size of Penalty (number of covariates)
  */
  LogisticLassoPenalty::LogisticLassoPenalty(STK::Real lambda)
                            : LassoPenalty(lambda)
                            , z_()
                            , p_y_(0)
                            , p_currentData_(0)
  {
  }


  /*
   * Copy constructor
   * @param penalty LassoPenalty object to copy
   */
  LogisticLassoPenalty::LogisticLassoPenalty(LogisticLassoPenalty const& penalty)
              : LassoPenalty(penalty)
              , z_(penalty.z())
              , p_y_(penalty.p_y())
              , p_currentData_(penalty.p_currentData())
  {
  }

  /*clone*/
  LogisticLassoPenalty* LogisticLassoPenalty::clone() const
  {
    return new LogisticLassoPenalty(*this);
  }

  /*
   * update the lasso penalty (for fixed sigma)
   * @param beta current estimates
   */
  void LogisticLassoPenalty::update(STK::CVectorX const& beta)
  {
    updatePenalty(beta);
    updateZ(beta);
  }


  void LogisticLassoPenalty::updateZ(STK::CVectorX const& beta)
  {

    STK::Law::Normal normal;

    for(int i = 1; i <= z_.sizeRows(); i++)
    {
      z_[i] = beta.dot(p_currentData_->row(i));

      if((*p_y_)[i] == 1)
      {
        z_[i] += normal.pdf(z_[i])/(1.-normal.cdf(-z_[i]));
      }
      else
      {
        z_[i] -= normal.pdf(z_[i])/normal.cdf(-z_[i]);
      }
    }

  }

}
