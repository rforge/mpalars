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
          , STK::Real eps =STK::Arithmetic<STK::Real>::epsilon(), bool verbose = false);

     //getters
      /**@return path of the coefficients*/
      inline Path const& path() const {return path_;}

      /**
       * @param i step
       * @return the Pathstate object : the state of the path at the step i
       */
      inline PathState const& path(int i) const {return path_.states(i);}

      /**
       * @param i index of the step
       * @param j index of the coefficients
       * @return the value of the j-th coefficient at the step i
       */
      inline STK::Real coefficient(int i,int j) const {return path_.varCoeff(i,j);}

      /**
       * @param i index of the step
       * @param j index of the coefficients
       * @return the value of the j-th coefficient at the step i
       */
      inline int varIdx(int i,int j) const {return path_.varIdx(i,j);}

      /**
       * @param i index of the step
       * @return the value of lambda at the i-th step
       */
      inline STK::Real lambda(int i) const {return path_.lambda(i);}

      /**
       * @return the vector of lambda
       */
      inline STK::Array2DVector<STK::Real> const lambda() const {return path_.lambda();}

      /** @return the historic of add and drop variable*/
      inline std::vector< std::pair<int,int> > evolution() const {return path_.evolution();}

      /**@return Number of step of the algorithm*/
      inline int step() const {return step_;}

      /** @return the intercept of the solution*/
      inline STK::Real mu() const {return mu_;}

      /** @return the ignored variable*/
      inline std::vector<bool> toIgnore() const {return toIgnore_;}



      //methods
      /** run lars algorithm*/
      void run();

      /**
       * predict the path for a ratio fraction = lambda/lambamax
       * @param X new data for predict the response
       * @param fraction real between 0 and 1 .
       * @return predicted response
       */
      void predict(STK::CArrayXX const& X, STK::Real fraction, STK::CVectorX &yPred);

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
      void computeAddSet(STK::Real Cmax, std::vector<int>& newId) const;

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
      void computeGi1(STK::Array2DVector<STK::Real> &Gi1, STK::Array2DVector<int> const& signC) const;

      /**
       * compute Cmax
       * @return Cmax the correlation max
       */
      STK::Real computeCmax();

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
       * @return gammaHat a STK::Real
       */
      STK::Real computeGamHat(STK::Real const& Aa, STK::CVectorX const& a, STK::Real Cmax) const;

      /**
       * Compute gammaTilde for the update of coefficient in drop case
       * @param w Aa*Gi1 @see computeGi1
       * @param idxMin we stock the index (in the activeVariable vector) of the variable with the min value
       * @return gammatilde a STK::Real
       */
      STK::Real computeGamTilde(STK::Array2DVector<STK::Real> const& w, int &idxMin) const;


      /**
       * Update the coefficient of the path
       * @param gamma gammaHat or gammaTilde
       * @param w Aa*Gi1 @see computeGi1
       * @param action a pair with first element is a bool (true for addcase, false for dropcase) and second the idx variable to drop/add
       * @param isAddCase true if we add a variable
       * @param dropId id top potentially drop
       */
      void updateBeta(STK::Real gamma, STK::Array2DVector<STK::Real> const& w, std::pair<bool,int> action, bool isAddCase, int dropId);


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
      bool firstStep(STK::Real &Cmax, std::vector<int> &newId, STK::Array2DVector<int> &signC, std::pair<bool,int> &action, STK::Real &Aa,
                     STK::Array2DVector<STK::Real> &Gi1, STK::Array2DVector<STK::Real> &w, STK::CVectorX &u, STK::CVectorX &a, STK::Real &gam);

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

      STK::Array2DVector<STK::Real> activeVariables_;
      int step_;
      STK::Real mu_;//B0
      STK::Real eps_;//eps for zero approximation
      bool verbose_;

      STK::Qr qrX_;//qr decomposition of Xi
      STK::Array2D<STK::Real> Xi_;//X for covariates of active Set
      STK::CVectorX c_;//correlation size p*1
  };

}//end namespace

#endif /* LARS_H_ */
