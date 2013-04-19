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
 * created on: 6 févr. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file functions.cpp
 *  @brief Contains utilities functions for lars algorithm.
 **/

#include "functions.h"

using namespace std;
using namespace MPA;

/*
 * (x1,y1) and (x2,y2) are points of an affine function, we wants to calculate the image of x3.
 * @param x1 abscissa of the first point
 * @param x2 abscissa of the second point
 * @param x3 abscissa of the point we want ordinate
 * @param y1 ordinate of the first point
 * @param y2 ordinate of the second point
 * @return ordinate of x3
 */
Real computeOrdinate(Real x1,Real x2,Real x3,Real y1,Real y2)
{
  return x1 + (y2-y1) * ((x3-x1)/(x2-x1)) ;
}

/*
 * Compute the coefficients for a given value of lambda
 * Use with move()
 * @param state1 state of a lars step
 * @param state2 state of the next lars step
 * @param evolution difference between the 2 lars step
 * @param lambda abscissa to compute ordinates
 * @return value of coefficients for lambda
 */
STK::Array2DVector< pair<int,Real> > computeCoefficients(PathState const& state1,PathState const& state2,pair<int,int> const& evolution, Real const& lambda)
{
  STK::Array2DVector< pair<int,Real> > coeff(max(state1.coefficients().size(),state2.coefficients().size()));

  if(evolution.second==0)
  {//no drop variable
    int j(0);
    for(j=0; j<state1.coefficients().size(); j++)
      coeff[j+1]=make_pair(state1.varIdx(j+1),computeOrdinate(state1.lambda(), state2.lambda(), lambda,state1.varCoeff(j+1), state2.varCoeff(j+1)));

    //add variable case
    if(evolution.first!=0)
      coeff[j+1]=make_pair( evolution.second, computeOrdinate(state1.lambda(),state2.lambda(),lambda,0.,state2.varCoeff(j+1)));
  }
  else
  {
    //delete variable case
    int i(1);
    //while we don't meet the delete variable, variable has the same index in the two sets
    while(evolution.second!=state1.coefficients(i).first)
    {
      coeff[i]=make_pair(state1.varIdx(i),computeOrdinate(state1.lambda(),state2.lambda(),lambda,state1.varCoeff(i),state2.varCoeff(i)));
      i++;
    }
    //compute coefficient for the delete variable
    coeff[i]=make_pair(state1.varIdx(i),computeOrdinate(state1.lambda(),state2.lambda(),lambda,state1.varCoeff(i),0.));
    i++;

    //compute coefficient for the other variable
    while(i<state1.coefficients().size()+1)
    {
      coeff[i]=make_pair(state1.varIdx(i),computeOrdinate(state1.lambda(),state2.lambda(),lambda,state1.varCoeff(i),state2.varCoeff(i-1)));
      i++;
    }

    //drop with an add variable
    if(evolution.first!=0)
      coeff[i]=make_pair(evolution.first, computeOrdinate(state1.lambda(), state2.lambda(), lambda, 0., state2.varCoeff(i)));

  }

  return coeff;

  return coeff;
}

void print(STK::Array2DVector< pair<int,Real> > const& state)
{
  for(int j(1);j<=state.size();j++)
    cout<<state[j].first<<"        ";
  cout<<endl;
  for(int j(1);j<=state.size();j++)
    cout<<state[j].second<<" ";
  cout<<endl;
}


void import(string adressFichier,int n,int p,STK::CArrayXX &data)
{
  ifstream flux(adressFichier.c_str());

  Real real;
  int i(1),j(1);

  if (flux)//si le fichier est ouvert
  {
     while(flux>>real)//on lit le fichier entier par entier
     {
       data(i,j)=real;
       j++;
       if(j>p)//bout de ligne, on ajoute les données et on passe à la ligne
       {
         j=1;
         i++;
       }
     }
  cout<<"Importation finie."<<endl<<endl;
  }
  else
    cout<<"Problème d'importation!"<<endl<<endl;
}

void import(string adressFichier,int n,STK::CVectorX &data)
{
  ifstream flux(adressFichier.c_str());

  Real real;
  int i(1);

  if (flux)//si le fichier est ouvert
  {
     while(flux>>real)//on lit le fichier entier par entier
     {
       data[i]=real;
       i++;
     }
  cout<<"Importation finie."<<endl<<endl;
  }
  else
    cout<<"Problème d'importation!"<<endl<<endl;
}
