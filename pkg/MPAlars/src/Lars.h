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
 * created on: 12 févr. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file Lars.h
 *  @brief we define the class @c Lars.
 **/


#ifndef LARS_H_
#define LARS_H_

#include "Path.h"


namespace MPA
{
  class Lars
  {
    public:
      //constructors
      /**
       * Constructor
       * @param X matrix of data, a row=a individual
       * @param y response
       */
      Lars(STK::CArrayXX const& X, STK::CVectorX const& y);

      /**
       * Constructor
       * @param X matrix of data, a row=a individual
       * @param y response
       * @param maxSteps number of maximum step to do
       * @param eps epsilon (for 0)
       * @param verbose if TRUE print some details
       */
      Lars( STK::CArrayXX const& X, STK::CVectorX const& y, int maxSteps
          , Real eps =STK::Arithmetic<Real>::epsilon(), bool verbose = false);

      //getters
      inline Path const& path() const {return path_;}
      inline PathState const& path(int i) const {return path_.states(i);}
      inline Real coefficient(int i,int j) const {return path_.varCoeff(i,j);}
      inline int varIdx(int i,int j) const {return path_.varIdx(i,j);}
      inline Real lambda(int i) const {return path_.lambda(i);}
      inline std::vector< std::pair<int,int> > evolution() const {return path_.evolution();}
      inline int step() const {return step_;}
      inline Real mu() const {return mu_;}
      inline std::vector<bool> toIgnore() const {return toIgnore_;}


      //methods
      /** run lars algorithm*/
      void run();

   protected:
      /**
       * initialization of algorithm
       */
      void initialization();

      /**
       * search non active variable with the greatest correlation
       * @param Cmax correlation max
       * @param newId a vector containing the index of variable to potentially add
       */
      void computeAddSet(Real Cmax, std::vector<int>& newId) const;

      /**
       * update the QR decomposition of Xi
       * @param idxVar index of active variable to add
       * @param signC sign of correlation of active variable
       * @param action a pair with first element is a bool (true for addcase, false for dropcase) and second the idx variable to drop/add
       */
      void updateR(int idxVar, STK::Array2DVector<int> &signC, std::pair<bool,int> &action);

      /**
       * compute inv(Xi'*Xi)*1 from qr decomposition
       * @param Gi1 for stock inv(Xi'*Xi)*1
       * @param signC sign of correlation of active variable
       */
      void computeGi1(STK::Array2DVector<Real> &Gi1, STK::Array2DVector<int> const& signC) const;

      /**
       * compute Cmax
       * @return Cmax the correlation max
       */
      Real computeCmax();

      /** dropStep
       * downdate qr decomposition,  X and signC
       * @param idxVar index of active variable to drop
       * @param signC sign of correlation of active variable
       */
      void dropStep(int idxVar, STK::Array2DVector<int> &signC);

      /**
       * Compute gammahat for the update of coefficient in add case
       * @param Aa norm of the inverse of G
       * @param a X' * equiangular vector
       * @param Cmax correlation max
       * @return gammaHat a real
       */
      Real computeGamHat(Real const& Aa, STK::CVectorX const& a, Real Cmax) const;

      /**
       * Compute gammaTilde for the update of coefficient in drop case
       * @param w Aa*Gi1 @see computeGi1
       * @param idxMin we stock the index (in the activeVariable vector) of the variable with the min value
       * @return gammatilde a real
       */
      Real computeGamTilde(STK::Array2DVector<Real> const& w, int &idxMin) const;


      /**
       * Update the coefficient of the path
       * @param gamma gammaHat or gammaTilde
       * @param w Aa*Gi1 @see computeGi1
       * @param action a pair with first element is a bool (true for addcase, false for dropcase) and second the idx variable to drop/add
       * @param isAddCase true if we add a variable
       * @param dropId id top potentially drop
       */
      void updateBeta(Real gamma, STK::Array2DVector<Real> const& w, std::pair<bool,int> action, bool isAddCase, int dropId);


      /**
       * first step
       * @param Cmax correlation max
       * @param newId vector of index of active variable to add
       * @param signC sign of correlation of active variable
       * @param action a pair with first element is a bool (true for addcase, false for dropcase) and second the idx variable to drop/add
       * @param Aa norm of the inverse of G
       * @param Gi1 for stock inv(Xi'*Xi)*1
       * @param w Aa*Gi1
       * @param u unit vector making equal angles with the column of Xi
       * @param a X' * equiangular vector       * @param Gi1
       * @param gam the step for update coefficients
       * @return
       */
      bool firstStep(Real &Cmax, std::vector<int> &newId, STK::Array2DVector<int> &signC, std::pair<bool,int> &action, Real &Aa,
                     STK::Array2DVector<Real> &Gi1, STK::Array2DVector<Real> &w, STK::CVectorX &u, STK::CVectorX &a, Real &gam);

      /**
       * updateR only for the first step
       * @see updateR
       */
      void firstUpdateR(int idxVar, STK::Array2DVector<int> &signC, std::pair<bool,int> &action);

    private:
      int n_;//number of individuals
      int p_;//number of variables
      int maxSteps_; // maximal number of steps

      STK::CArrayXX X_;//covariate size n*p
      STK::CVectorX y_;//response size p*1
      STK::CVectorX muX_;//mean of each covariate of X

      Path path_;

      std::vector<bool> isActive_;
      std::vector<bool> toIgnore_;//id to ignore because it causes singularity

      int nbActiveVariable_;
      int nbIgnoreVariable_;

      STK::Array2DVector<Real> activeVariables_;
      int step_;
      Real mu_;//B0
      Real eps_;//eps for zero approximation
      bool verbose_;

      STK::Qr qrX_;//qr decomposition of Xi
      STK::Array2D<Real> Xi_;//X for covariates of active Set
      STK::CVectorX c_;//correlation size p*1
  };

}//end namespace

#endif /* LARS_H_ */
