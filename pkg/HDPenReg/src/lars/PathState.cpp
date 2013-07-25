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

    Contact : Serge.Iovleff@stkpp.org
*/

/*
 * Project:  MPAGenomics::
 * created on: 5 févr. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file PathState.cpp
 *  @brief Code of methods associates to @c PathState.
 **/

#include "PathState.h"

using namespace std;
using namespace STK;

namespace HD
{
  //constructors
  /* default constructor*/
  PathState::PathState() : l1norm_(0.)
  {}

  /*
   * constructor with reserve size
   * @param nbMaxVariable maximal number of variable to potentially stock
   */
  PathState::PathState(int nbMaxVariable) : l1norm_(0.)
  {coefficients_.reserveCols(nbMaxVariable);}

  /* update coefficients with values specified in parameters
   * @param indexVariables vector containing the index of active variables
   * @param coefficients vector containing the value of estimates for active variables
   */
  void PathState::update(STK::Array2DVector<int> const& indexVariables,STK::Array2DVector<Real> const& coefficients)
  {
    //resize the container
    coefficients_.resize(coefficients.size());
    l1norm_=0;
    //fill the container and compute l1norm
    for(int i=1; i<=coefficients.size(); i++)//array index starts at 1
    {
      coefficients_[i]=make_pair(indexVariables[i],coefficients[i]);
      l1norm_ += std::abs(coefficients[i]);
    }
  }

  /*print coefficients*/
  void PathState::printCoeff() const
  {
    for(int j(1);j<=coefficients_.size();j++)
      cout<<coefficients_[j].first<<"        ";
    cout<<endl;
    for(int j(1);j<=coefficients_.size();j++)
      cout<< coefficients_[j].second<<" ";
    cout<<endl;
  }

  /*
   * update of the coefficients of the previous step
   * @param w direction of the update
   * @param gamma step of the update
   */
  void PathState::update(STK::Array2DVector<Real> const& w, Real gamma)
  {
    l1norm_=0;
    for(int i=1; i<=(*this).size(); i++)
    {
      coefficients_[i].second+=gamma*w[i];
      l1norm_ += std::abs(coefficients_[i].second);
    }
  }

  /*
   * update of the coefficients of the previous state with a variable to drop and a variable to add
   * @param w direction of the update
   * @param gamma step of the update
   * @param addIdxVar index of the variable to add
   * @param dropIdx index (in the vector of coefficients of the previous step) of the variable to delete
   */
  void PathState::addWithDropUpdate(STK::Array2DVector<Real> const& w, Real gamma, int addIdxVar, int dropIdx)
  {
    //update the other variable
    l1norm_=0;
    if(dropIdx!=1)
      for(int i=1; i<dropIdx; i++)
      {
        coefficients_[i].second += gamma*w[i];
        l1norm_ += std::abs(coefficients_[i].second);
      }

    if(dropIdx!=coefficients_.size())
      for(int i=dropIdx+1; i<=coefficients_.size(); i++)
      {
        coefficients_[i].second += gamma*w[i];
        l1norm_ += std::abs(coefficients_[i].second);
      }

    //add the new variable
    coefficients_.pushBack();
    coefficients_.back()=make_pair(addIdxVar,gamma*w.back());
    l1norm_ += std::abs(coefficients_.back().second);

    //delete the variable to delete
    coefficients_.erase(dropIdx,1);
  }

  /*
   * update of the coefficients of the previous state with a variable to drop
   * @param w direction of the update
   * @param gamma step of the update
   * @param dropIdx index (in the vector of coefficients of the previous step) of the variable to delete
   */
  void PathState::dropAfterDropUpdate(STK::Array2DVector<Real> const& w, Real gamma, int dropIdx)
  {
    l1norm_=0;
    //update the other coefficient
    if(dropIdx!=1)
      for(int i=1; i<dropIdx; i++)
      {
        coefficients_[i].second += gamma*w[i];
        l1norm_ += std::abs(coefficients_[i].second);
      }

    if(dropIdx!=coefficients_.size())
      for(int i=dropIdx+1; i<=coefficients_.size(); i++)
      {
        coefficients_[i].second += gamma*w[i];
        l1norm_ += std::abs(coefficients_[i].second);
      }

    //delete the variable to delete
    coefficients_.erase(dropIdx,1);
  }

  /*
   * update of the coefficients of the previous state with a new variable
   * @param w direction of the update
   * @param gamma step of the update
   * @param addIdxVar index of the variable to add
   */
  void PathState::addUpdate(STK::Array2DVector<Real> const& w, Real gamma, int addIdxVar)
  {
    //update previous coefficients
    l1norm_=0;
    for(int i=1; i<=coefficients_.size(); i++)
    {
      coefficients_[i].second += gamma * w[i];
      l1norm_ += std::abs(coefficients_[i].second);
    }

    coefficients_.pushBack(1);
    coefficients_.back()=make_pair(addIdxVar,gamma * w.back());
    l1norm_ += std::abs(coefficients_.back().second);
  }

}//end namespace



