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
 * Project:  stkpp::tests
 * created on: 8 août 2011
 * Purpose:  test the Normal and MultiNormal classes.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file testJointModels.cpp
 *  @brief In this file we test the joint Statistical models.
 **/

#include "../../include/STKpp.h"

using namespace STK;

// initialize static member generator
RandBase Law::ILawBase::generator;

/* main. */
int main(int argc, char *argv[])
{
  int N = (argc < 2) ? 200 : int(atoi(argv[1]));

  stk_cout << _T("\n\n");
  stk_cout << _T("+++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  stk_cout << _T("+ Test JointBernoulliModel                          +\n");
  stk_cout << _T("+++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  stk_cout << _T("\n\n");

  // generate Bernoulli joint law
  Law::JointBernoulli<Array2DPoint<Binary> > bLaw(20);
  for (Integer j= bLaw.JointLaw().firstIdx(); j <= bLaw.JointLaw().lastIdx(); ++j)
  { bLaw.setProb(j, Law::ILawBase::generator.randUnif());}
  // generate data set
  Array2D<Binary> bData(N, 20);
  for (Integer i= bData.firstIdxRows(); i <= bData.lastIdxRows(); ++i)
  {
    Array2DPoint<Binary> row(bData.row(i), true);
    bLaw.rand(row);
  }
  // run model
  JointBernoulliModel<Array2D<Binary> > bModel(bData);
  bModel.run();
  for (Integer j= bLaw.JointLaw().firstIdx(); j <= bLaw.JointLaw().lastIdx(); ++j)
  {
    stk_cout << _T("j= ") << j << _T(". True parameter= ") << bLaw.prob(j);
    stk_cout << _T(". Estimated parameter= ") << bModel.p_param()->prob(j) << _T("\n");
  }

  stk_cout << _T("\n\n");
  stk_cout << _T("+++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  stk_cout << _T("+ Successful completion of testing for              +\n");
  stk_cout << _T("+ JointBernoulliModel.                              +\n");
  stk_cout << _T("+++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  stk_cout << _T("\n\n");

  stk_cout << _T("\n\n");
  stk_cout << _T("+++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  stk_cout << _T("+ Test JointGaussianModel                           +\n");
  stk_cout << _T("+++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  stk_cout << _T("\n\n");

  // run model
  // generate Bernoulli joint law
  Law::JointNormal<Array2DPoint<Real> > nLaw(20);
  for (Integer j= nLaw.JointLaw().firstIdx(); j <= nLaw.JointLaw().lastIdx(); ++j)
  {
    nLaw.setMu(j, Law::ILawBase::generator.randGauss());
    nLaw.setSigma(j, Law::ILawBase::generator.randExp());
  }
  // generate data set
  Array2D<Real> nData(N, 20);
  for (Integer i= nData.firstIdxRows(); i <= nData.lastIdxRows(); ++i)
  {
    Array2DPoint<Real> row(nData.row(i), true);
    nLaw.rand(row);
  }
  // run model
  JointGaussianModel<Array2D<Real> > nModel(nData);
  nModel.run();
  for (Integer j= nLaw.JointLaw().firstIdx(); j <= nLaw.JointLaw().lastIdx(); ++j)
  {
    stk_cout << _T("j= ") << j << _T(". True mu= ") << nLaw.mu(j);
    stk_cout << _T(". Estimated mu= ") << nModel.p_param()->mu(j) << _T("\n");
    stk_cout << _T("j= ") << j << _T(". True sigma= ") << nLaw.sigma(j);
    stk_cout << _T(". Estimated sigma= ") << nModel.p_param()->sigma(j) << _T("\n\n");
  }

  stk_cout << _T("\n\n");
  stk_cout << _T("+++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  stk_cout << _T("+ Successful completion of testing for              +\n");
  stk_cout << _T("+ JointGaussian Model.                              +\n");
  stk_cout << _T("+++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  stk_cout << _T("\n\n");

  return 0;
}

