/*--------------------------------------------------------------------*/
/*     Copyright (C) 2013-2013  Serge Iovleff, Quentin Grimonprez

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

    Contact : quentin.grimonprez@inria.fr
*/

/*
 * Project:  MPAGenomics::
 * created on: 3 avr. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file Fusion.cpp
 *  @brief In this file, we define the methods of the @c Fusion class .
 **/

#include "Fusion.h"

using namespace std;
using namespace STK;


namespace HD
{
  //constructors
  
  /*
   * constructor
   * @param X matrix of size n*p, a row contains the values of each covariate for an individual.
   * @param y vector of length n containing the response
   */
  Fusion::Fusion(STK::CArrayXX const& X, STK::CVectorX const& y, bool intercept)
                : X_(X),
                  y_(y),
                  eps_(STK::Arithmetic<Real>::epsilon()),
                  path_(maxSteps_),
                  toIgnore_(),
                  intercept_(intercept)
  {
    maxSteps_ = 3*min(X.sizeRows(),X.sizeCols());
    computeZ();
  }
  
  /*
   *
   * @param X matrix of size n*p, a row contains the values of each covariate for an individual.
   * @param y vector of length n containing the response
   * @param maxSteps number of maximum step to do
   * @param eps epsilon (for 0)
   */
  Fusion::Fusion( STK::CArrayXX const& X, STK::CVectorX const& y, int maxSteps, bool intercept, Real eps)
                : X_(X),
                  y_(y),
                  maxSteps_(maxSteps),
                  eps_(eps),
                  path_(maxSteps),
                  toIgnore_(),
                  intercept_(intercept)
  {
    computeZ();
  }
  
  //destructors
  
  
  //methods
  /*
   * change X in Z=X*L^-1 (L^-1 = lower triangular matrix of 1)
   */
  void Fusion::computeZ()
  {
    int p=X_.sizeCols(), n=X_.sizeRows();
  
    for(int i=p-1; i>=1; i--)
    {
      for(int j=1; j<=n; j++ )
        X_(j,i) += X_(j,i+1);
  
    }
  //    X_.col(i) += X_.col(i-1);
  
  }
  
  /*
   * run the lars algorithm for solving the fusion problem on Z=X*L^-1 (L^-1 = lower triangular matrix of 1)
   */
  void Fusion::run()
  {
    //run lars algorithm on Z
    Lars lars(X_, y_, maxSteps_, intercept_, eps_);
    lars.run();
  
    //get the solution path
    path_=lars.path();
    step_=lars.step();
    mu_=lars.mu();
    muX_=lars.muX();
    toIgnore_=lars.toIgnore();
    msg_error_=lars.msg_error();
  
    //transform the coefficient B in L^-1 * B
    //int p=X_.sizeCols();
    //if(step_>=2)
    //  path_.transform2fusion(p);
  }

}//end namespace
