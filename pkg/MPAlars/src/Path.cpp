/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2013  Serge Iovleff

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
 * created on: 4 févr. 2013
 * Author:   grimonprez, serge.iovleff@stkpp.org
 **/

/** @file Path.cpp
 *  @brief In this file .
 **/

#include "functions.h"
#include "Path.h"
#include <cmath>
#include <iostream>

using namespace std;

namespace MPA
{
  //Constructors
  Path::Path(int maxSizePath)
  {
    states_.reserve(maxSizePath);
    evolution_.reserve(maxSizePath);
    states_.push_back(PathState());
  }

  //Destructors

  //Methods

  /* @brief Add coefficients of a LARS step to the actual path
   * @param indexVariables Array2DVector containing the index of active variables
   * @param coefficients Array2DVector containing the value of estimates for the active variables
   * @param idxVarAdd index of the new variable (0 if no new variable)
   * @param idxVarDrop index of the delete variable variable (0 if no delete variable)
   */
  void Path::addCoeff(STK::Array2DVector<int> const& indexVariables,STK::Array2DVector<Real> const& coefficients,int idxVarAdd,int idxVarDrop)
  {
    states_.push_back(PathState());
    states_.back().update(indexVariables,coefficients);
    evolution_.push_back(make_pair(idxVarAdd,idxVarDrop));
  }


  /* @brief get coefficient associates to a lambda value.
   *  @param lambda is the norm value for which we want values of coefficient
   *  @return a vector containing pair<int,double>=(index of non zero coefficient,coefficient)
   */
  STK::Array2DVector< pair<int,Real> > Path::coeff(Real lambda) const
  {
    STK::Array2DVector< pair<int,Real> > coeff;
    if(lambda!=0)
    {
      //search off the interval of m_lambda containing lambda
      int indexLambda=0;
      while( ( states_[indexLambda].lambda() < lambda ) && ( indexLambda < states_.size()-1 ) )
        indexLambda++;

      if(lambda==states_[indexLambda].lambda())//if lambda is equal to an actual lambda
        coeff=states_[indexLambda].coefficients();
      else
      {
        if( indexLambda == states_.size()-1 )//last m_lambda, we return the last coefficient
          coeff=states_[indexLambda].coefficients();
        else
          coeff=computeCoefficients(states_[indexLambda-1],states_[indexLambda],evolution_[indexLambda-1],lambda);
      }
    }
    return coeff;
  }

  /*
   * update of the coefficients of the previous state with a new variable
   * @param w direction of the update
   * @param gamma step of the update
   * @param addIdxVar index of the variable to add
   */
  void Path::addCaseUpdate(Real gamma, STK::Array2DVector<Real> const& w, int addIdxVar)
  {
    //copy the previous states
    states_.push_back(states_.back());

    //update of evolution with the new index
    evolution_.push_back(make_pair(addIdxVar,0));

    //update of the coefficients
    states_.back().addUpdate(w,gamma,addIdxVar);
  }

  /*
   * update of the coefficients of the previous step
   * @param w direction of the update
   * @param gamma step of the update
   */
  void Path::update(Real gamma, STK::Array2DVector<Real> const& w)
  {
    //copy the previous states
    states_.push_back(states_.back());

    //update of evolution with the new index
    evolution_.push_back(make_pair(0,0));

    //update of the coefficients
    states_.back().update(w,gamma);
  }

  /*
   * update of the coefficients of the previous state with a variable to drop and a variable to add
   * @param w direction of the update
   * @param gamma step of the update
   * @param addIdxVar index of the variable to add
   * @param dropIdxVar index of the delete variable
   * @param dropIdx index (in the vector of coefficients of the previous step) of the variable to delete
   */
  void Path::addWithDropCaseUpdate(Real gamma, STK::Array2DVector<Real> const& w, int addIdxVar, int dropIdxVar, int dropIdx)
  {
    //copy the previous states
    states_.push_back(states_.back());

    //update of evolution with the new index
    evolution_.push_back(make_pair(addIdxVar,dropIdxVar));

    //update of the coefficients
    states_.back().addWithDropUpdate(w, gamma, addIdxVar, dropIdx);
  }

  /*
   * update of the coefficients of the previous state with a variable to drop
   * @param w direction of the update
   * @param gamma step of the update
   * @param dropIdxVar index of the delete variable
   * @param dropIdx index (in the vector of coefficients of the previous step) of the variable to delete
   */
  void Path::dropAfterDropCaseUpdate(Real gamma, STK::Array2DVector<Real> const& w, int dropIdxVar, int dropIdx)
  {
    //copy the previous states
    states_.push_back(states_.back());

    //update of evolution with the new index
    evolution_.push_back(make_pair(0,dropIdxVar));

    //update of the coefficients
    states_.back().dropAfterDropUpdate(w,gamma,dropIdx);
  }

  /*
   * transform the path from lars problem to fusion problem
   * @param p Number of variables
   */
  void Path::transform2fusion(int const& p)
  {
    STK::Array2DVector<int> order(1,1);// the first element contains the index (in the vector
    order.reserveCols(states_.size());
    STK::Array2DVector< pair<int,Real> > coeffTemp;
    coeffTemp.reserveCols(states_.size());

    Real lambda;
    int drop(1);

    //step=1. only 1 variable, so we only update lambda
    lambda = abs((p-states_[1].varIdx(1))*states_[1].varCoeff(1));
    states_[1].setLambda(lambda);

    for(int step=2; step<states_.size(); step++)
    {
      //if a variable is drop at this step
      if(evolution_[step-1].second!=0)
      {
        //we delete the position found
        for(int j=1; j<=order.size(); j++)
        {
          if( order[j]==drop )
          {
            order.erase(j,1);
            break;
          }
        }

        for(int j=1; j<=order.size(); j++)
        {
          if( order[j]>drop )
            order[j]--;
        }
      }

      //if a variable is add at this step
      if(evolution_[step-1].first!=0)
      {
        //update order
        for(int j=1; j<states_[step].size(); j++)
        {
          //we look after the first index greater than the new index
          if( states_[step].varIdx(states_[step].size()) < states_[step].varIdx(order[j]) )
          {
            order.insertElt(j,1);
            order[j]=states_[step].size();

            break;
          }
        }
        if(order.size()<states_[step].size())
        {
          order.pushBack(1);
          order.back()=states_[step].size();
        }
      }

      //if there is a drop variable at the next step, we search its position in the actual coefficients before coefficients are sorting
      if(step<states_.size()-1)
      {
        if(evolution_[step].second!=0)
        {
          for(int i=1; i<=states_[step].size(); i++)
          {
            if(states_[step].varIdx(i)==evolution_[step].second)
            {
              drop=i;
              break;
            }
          }
        }
      }

      //transform to fusion coefficient
      coeffTemp=states_[step].coefficients();
      lambda=0;

      coeffTemp[1]=states_[step].coefficients(order[1]);

      for(int i=2; i<states_[step].size(); i++)
      {
        //we sort the index
        coeffTemp[i]=states_[step].coefficients(order[i]);

        //compute coefficients of fusion
        coeffTemp[i].second += coeffTemp[i-1].second;

        //update lambda
        lambda += abs((coeffTemp[i].first-coeffTemp[i-1].first)*coeffTemp[i-1].second);
      }

      coeffTemp.back()=states_[step].coefficients(order.back());
      coeffTemp.back().second += coeffTemp[coeffTemp.size()-1].second;
      lambda += abs((p-coeffTemp.back().first)*coeffTemp.back().second);

      states_[step].setCoefficients(coeffTemp);
      states_[step].setLambda(lambda);

    }
  }
}//end namespace




