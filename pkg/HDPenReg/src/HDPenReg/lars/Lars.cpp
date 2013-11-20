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
 * created on: 13 févr. 2013
 * Author:   Quentin Grimonprez
 **/

/** @file Lars.cpp
 *  @brief In this file, methods associates to @c Lars.
 **/

#include <cmath>
#include <algorithm>
#include <limits>
#include "Lars.h"
#include "functions.h"

using namespace STK;
using namespace std;

namespace HD
{
  //Constructors
  /*
    * Constructor
    * @param X matrix of data, a row=a individual
    * @param y response
    */
  Lars::Lars( CArrayXX const& X,CVectorX const& y, bool intercept)
            : n_(X.rows().size())
            , p_(X.cols().size())
            , maxSteps_(3*min(n_,p_))
            , X_(X), y_(y)
            , muX_(p_)
            , path_(maxSteps_)
            , isActive_(p_, false)
            , toIgnore_(p_,false)
            , nbActiveVariable_(0)
            , nbIgnoreVariable_(0)
            , activeVariables_(0)
            , step_(0)
            , mu_()
            , eps_(STK::Arithmetic<Real>::epsilon())
            , qrX_()
            , Xi_(n_,1)
            , c_()
            , intercept_(intercept)
            , msg_error_()
  { initialization();}

  /*
   *
   * @param X matrix of data, a row=a individual
   * @param y response
   * @param maxStep number of maximum step to do
   * @param eps epsilon (for 0)
   */
  Lars::Lars( CArrayXX const& X,CVectorX const& y, int maxSteps, bool intercept, Real eps)
            : n_(X.rows().size())
            , p_(X.cols().size())
            , maxSteps_(maxSteps)
            , X_(X)
            , y_(y)
            , muX_(p_)
            , path_(maxSteps_)
            , isActive_(p_, false)
            , toIgnore_(p_,false)
            , nbActiveVariable_(0)
            , nbIgnoreVariable_(0)
            , activeVariables_(0)
            , step_(0)
            , mu_()
            , eps_(eps)
            , qrX_()
            , Xi_(n_,1)
            , c_()
            , intercept_(intercept)
            , msg_error_()
  { initialization();}

  //Methods
  /* initialization of algorithm
   */
  void Lars::initialization()
  {
    if(intercept_)
    {
      //we center y
      mu_=y_.sum()/n_;
      y_-=mu_;

      for(int j=1;j <= p_; j++)
      {
        muX_[j] = X_.col(j).sum()/n_;
        X_.col(j) -= muX_[j];
      }
    }
    else
    {
      mu_=0;
      muX_.zeros();
    }


    c_=X_.transpose()*y_;
    Xi_.reserveCols(min(n_,p_));
  }


  /*
   * compute Cmax the maximum correlation
   */
  Real Lars::computeCmax()
  {
    Real Cmax = 0.;
    for(int i = 1; i <= p_; i++)
    {
      if(!isActive_[i-1])
        Cmax = max(Cmax, std::abs(c_[i]));
    }
    return Cmax;
  }

  /*
   * search non active variable with the greatest correlation
   * @param Cmax correlation max
   * @param newId a vector containing the index of variable to potentially add
   */
  void Lars::computeAddSet(Real Cmax, vector<int>& newId ) const
  {
    for(int i=0; i<p_; i++)
    {
      if(!isActive_[i])//update with only variable non active
      {

//        stk_cout<<(abs(c_[i+1])-Cmax)<<"  ";

        //don't care about ignore variable
        if( !toIgnore_[i] & (abs(c_[i+1]) >= Cmax-eps_) )
        {
          newId.push_back(i+1);
        }
      }
    }
//    stk_cout<<std::endl;
  }


  /*
   * update the QR decomposition of Xi
   * @param idxVar index of active variable to add
   * @param signC sign of correlation of active variable
   * @param action a pair with first element is a bool (true for addcase, false for dropcase) and second the idx variable to drop/add
   */
  void Lars::updateR(int idxVar,Array2DVector<int> &signC, pair<bool,vector<int> > &action)
  {
    //update Xi_
    Xi_.pushBackCols(1);
    Xi_.col(Xi_.lastIdxCols()) = X_.col(idxVar);

    //update the QR decomposition
    qrX_.pushBackCol(Xi_.col(Xi_.lastIdxCols()));

#ifdef VERBOSE
      cout<<"Step "<<step_<<" : Variable "<< idxVar<<" added"<<endl;
//      cout<<std::abs(qrX_.R()( min(n_,nbActiveVariable_+1), nbActiveVariable_+1) )<<std::endl;
#endif

    //check if the variable added is not colinear with an other
    if( std::abs(qrX_.R()( min(n_,nbActiveVariable_+1), nbActiveVariable_+1) ) < eps_ )
    {
      //we cancel the add of the variable in the qr decomposition
      qrX_.popBackCols();

#ifdef VERBOSE
        cout<<"Step "<<step_<<" : Variable "<< idxVar<<" dropped (colinearity)"<<endl;
#endif


      toIgnore_[idxVar-1]=true;//the variable is add to the ignore set
      nbIgnoreVariable_++;
      Xi_.popBackCols(1);
    }
    else
    {
      //update action

      action.second.push_back(idxVar);
      action.first = true;

      //add the index to the set of active variable
      activeVariables_.pushBack(1);
      activeVariables_.back() = idxVar;

      nbActiveVariable_++;
      isActive_[idxVar-1] = true;

      //compute signC
      signC.pushBack(1);
      signC[nbActiveVariable_] = ( c_[idxVar] > 0 ) ? 1 : -1;
//      cout<<"add "<<idxVar<<"  cor "<<c_[idxVar]<<" sign "<<signC[nbActiveVariable_]<<endl;

    }

  }


  /*
   * compute inv(Xi'*Xi)*1 from qr decomposition
   * @param Gi1 for stock inv(Xi'*Xi)*1
   * @param signC sign of correlation of active variable
   */
  void Lars::computeGi1(Array2DVector<Real> &Gi1, Array2DVector<int> const& signC) const
  {
    CVectorX v(nbActiveVariable_);

    Gi1.resize(nbActiveVariable_);
    Gi1=0;

    for(int i=1; i<=nbActiveVariable_; i++)
    {
      //resolve R'*v=signC(i)*e_i
      v=0.;
      v[i]=signC[i];

      v[1] /= qrX_.R()(i,i);

      for(int j=2; j<=nbActiveVariable_; j++)
      {
        for(int k=1; k<=(j-1); k++)
          v[j] -= qrX_.R()(k,j)*v[k];
        v[j]/=qrX_.R()(j,j);
      }

      //resolve R*signC*z=v
      v[nbActiveVariable_] *= signC[nbActiveVariable_]/qrX_.R()(nbActiveVariable_,nbActiveVariable_);

      for(int j=nbActiveVariable_-1; j>0; j--)
      {
        for(int k=j+1; k <= nbActiveVariable_; k++)
          v[j] -= signC[k] * qrX_.R()(j,k) * v[k];
        v[j] /= signC[j] * qrX_.R()(j,j);
      }
      Gi1+=v;
    }
  }

  /*
   * Compute gammahat for the update of coefficient in add case
   * @param Aa norm of the inverse of G
   * @param a X' * equiangular vector
   * @param Cmax correlation max
   * @return gammaHat a real
   */
  Real Lars::computeGamHat(Real const& Aa, CVectorX const& a, Real Cmax) const
  {
    Real gamHat(Cmax/Aa), gam(0);
    for(int i=1; i<=p_; i++)
    {
      //only for the non active variable and non ignored variable
      if(!isActive_[i-1] && !toIgnore_[i-1])
      {
        //gamma is the min only on the positive value of gam1 and gam2
        if(Aa!=a[i])
        {
          if( (gam = (Cmax-c_[i]) / (Aa-a[i]) ) > eps_ ) gamHat=min(gamHat,gam);
        }

        if(Aa!=-a[i])
        {
          if( ( gam = (Cmax+c_[i]) / (Aa+a[i]) ) > eps_ ) gamHat=min(gamHat,gam);
        }
      }
    }
    return gamHat;
  }

  /*
   * Compute gammaTilde for the update of coefficient in drop case
   * @param w Aa*Gi1 @see computeGi1
   * @param idxMin we stock the index (in the activeVariable vector) of the variable with the min value
   * @return gammatilde a real
   */
  Real Lars::computeGamTilde(Array2DVector<Real> const& w,vector<int> &idxMin) const
  {
    Real gamTilde(std::numeric_limits<Real>::max()),gam(0);
    idxMin.erase(idxMin.begin(),idxMin.end());
    for(int i=1; i <= path_.lastState().size(); i++)
    {
      if(w[i]) gam = -path_.lastVarCoeff(i)/w[i];

      //we search the minimum only on positive value
      if(gam > eps_)
      {
        if(gam < gamTilde)
        {
          gamTilde = gam;
          idxMin.erase(idxMin.begin(),idxMin.end());
        }
        if(gam==gamTilde)
          idxMin.push_back(i);
      }
    }

    return gamTilde;
  }

  /*
   * Update the coefficient of the path
   * @param gamma gammaHat or gammaTilde
   * @param w Aa*Gi1 @see computeGi1
   * @param action a pair with first element is a bool (true for addcase, false for dropcase) and second the idx variable to drop/add
   * @param isAddCase true if we add a variable
   * @param dropId id top potentially drop
   */
  void Lars::updateBeta(Real gamma, STK::Array2DVector<Real> const& w, pair<bool,vector<int> > action, bool isAddCase, vector<int> dropId)
  {
    if( (action.first) && (isAddCase) )
    {
      //add situation
      path_.addCaseUpdate(gamma,w,action.second);
    }
    else
    {
      if( (action.first) && (!isAddCase) )
      {
        //add situation with a drop
        vector<int> drop(dropId.size());
        for(int i = 0; i < (int) dropId.size(); i++)
          drop[i] = activeVariables_[dropId[i]];

        path_.addWithDropCaseUpdate(gamma,w,action.second,drop,dropId);
      }
      else
      {
        if((!action.first) && (isAddCase))
        {
          //update after a drop situation
          path_.update(gamma,w);
        }
        else
        {
          vector<int> drop(dropId.size());
          for(int i = 0; i < (int) dropId.size(); i++)
            drop[i] = activeVariables_[dropId[i]];

          //drop step after a drop step
          path_.dropAfterDropCaseUpdate(gamma,w,drop,dropId);
        }
      }
    }
  }


  /* dropStep
   * downdate qr decomposition,  X and signC
   * @param idxVar index of active variable to drop
   * @param signC sign of correlation of active variable
   */
  void Lars::dropStep(vector<int> const& idxVar,Array2DVector<int> &signC)
  {
#ifdef VERBOSE
      for(int i = 0; i < (int) idxVar.size(); i++)
        cout<<"Step "<<step_+1<<" : Variable "<< activeVariables_[idxVar[i]]<<" dropped"<<endl;
#endif

      //idxVar are sort in the increasing order, we erase first the element with the greater index
    for(int i = idxVar.size()-1; i >= 0 ; i--)
    {
      //downdate R
      qrX_.eraseCol(idxVar[i]);

      //downdate Xi_
      Xi_.eraseCols(idxVar[i],1);
      signC.erase(idxVar[i],1);
      isActive_[activeVariables_[idxVar[i]]-1] = false;
      activeVariables_.erase(idxVar[i],1);
    }

  }

  /*
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
  bool Lars::firstStep(Real &Cmax, vector<int> &newId, Array2DVector<int> &signC, pair<bool, vector<int> > &action,
                       Real &Aa, Array2DVector<Real> &Gi1, Array2DVector<Real> &w, CVectorX &u, CVectorX &a, Real &gam)
  {
    step_++;

    //computation of correlation
    Cmax=computeCmax();
    if( Cmax < eps_*100)
    {
      step_--;
#ifdef VERBOSE
      std::cout << "Correlation max is equal to 0.";
#endif

      return false;
    }

    newId.resize(0);
    computeAddSet(Cmax, newId);

    if(newId.size()==0)
    {
      step_--;
#ifdef VERBOSE
      std::cout << "No variable selected for add in the add step."<<std::endl;
#endif

      return false;
    }

    addCmax(Cmax);

    for (vector<int>::iterator it = newId.begin() ; it != newId.end(); it++)
      firstUpdateR(*it,signC,action);

    //compute the inverse of G
    computeGi1(Gi1,signC);

    //compute Aa
    Aa = 1/sqrt(Gi1.sum());

    //compute w
    w = Gi1 * signC * Aa;

    //compute equiangular vector
    u = Xi_ * w;

    //computation of gamma hat
    //if the number of active variable is equal to the max number authorized, we don't search a new index
    if( nbActiveVariable_ == min(n_-1,p_-nbIgnoreVariable_) )
       gam=Cmax/Aa;
    else
    {
      //computation of a
      a = X_.transpose() * u;

      //computation of gamma hat
      gam = computeGamHat(Aa,a,Cmax);
    }

    //update beta
    vector<int> vide;
    updateBeta(gam,w,action,true,vide);

    //update of c
    c_ -= (X_.transpose() * u) * gam;

    return true;
  }

  /*
   * updateR only for the first step
   * @see updateR
   */
  void Lars::firstUpdateR(int idxVar,Array2DVector<int> &signC, pair<bool,vector<int> > &action)
  {
    //create Xi_ and qrXi_
    Xi_.col(Xi_.lastIdxCols()) = X_.col(idxVar);
    qrX_ = Qr(Xi_);//creation of the first QR decomposition

#ifdef VERBOSE
      cout<<"Step 1 : Variable "<< idxVar<<" added"<<endl;
#endif

    //check if the variable added is not colinear with an other
    if( std::abs(qrX_.R()( min(n_,nbActiveVariable_+1), nbActiveVariable_+1) ) < eps_ )
    {
      qrX_.popBackCols();
#ifdef VERBOSE
        cout<<"Step 1 : Variable "<< idxVar<<" dropped (colinearity)"<<endl;
#endif

      toIgnore_[idxVar-1] = true;//the variable is add to the ignore set
      nbIgnoreVariable_++;
      Xi_.popBackCols(1);
    }
    else
    {
      //Add the idx to the active set
      activeVariables_.pushBack(1);
      activeVariables_.back() = idxVar;

      nbActiveVariable_++;
      isActive_[idxVar-1] = true;

      //compute signC
      signC.pushBack(1);
      signC[nbActiveVariable_] = ( c_[idxVar] > 0 ) ? 1 : -1;

      action.first = true;
      action.second.push_back(idxVar);
    }

  }

  /** run lars algorithm*/
  void Lars::run()
  {
#ifdef VERBOSE
      cout<<"######################################"<<endl;
      cout<<"########### LARS ALGORITHM ###########"<<endl;
      cout<<"######################################"<<endl<<endl;
#endif

    Chrono::start();

    //initialization();

    bool isAddCase(true), continuer;
    vector<int> dropId;
    Real Aa(0),gam(0),gammaTilde(0),Cmax(0);

    Array2DVector<Real> Gi1(1), w;
    Array2DVector<int> signC;
    vector<int> newId;
    CVectorX a(p_,0),u(n_,0);
    pair<bool,vector<int> > action;

    newId.reserve(p_);
    Gi1.reserveCols(min(n_-1,p_));
    w.reserveCols(min(n_-1,p_));
    signC.reserveCols(min(n_-1,p_));

    continuer=firstStep(Cmax,newId,signC,action,Aa,Gi1,w,u,a,gam);
    //we stop, if we reach maxStep or if there is no more variable to add
    if(continuer)
    {
      STK::Real oldCmax;
      while( (step_< maxSteps_) && ( nbActiveVariable_ < min( n_-1, (p_-nbIgnoreVariable_) ) ) )
      {
        step_++;
        oldCmax = Cmax;
        //computation of correlation
        Cmax = computeCmax();

        if( Cmax < eps_*100)
        {
          step_--;
#ifdef VERBOSE
          std::cout << "Correlation max is equal to 0."<<std::endl;
#endif
          msg_error_ = "Correlation max is equal to 0.";
          break;
        }

        //if correlation max increased, we stop, Cmax must decreased
        if( Cmax > oldCmax)
        {
//          stk_cout<<Cmax<<endl;
          step_--;
#ifdef VERBOSE
          std::cout << "Correlation max has increased."<<std::endl;
#endif
          msg_error_ = "Correlation max has increased";
          break;
        }

        //add case : update of QR decomposition, active set and X'*X
        if(isAddCase)
        {
          newId.resize(0);
          computeAddSet(Cmax, newId);

          if(newId.size() == 0)
          {
            step_--;
#ifdef VERBOSE
            std::cout << "No variable selected for add in the add step."<<std::endl;
#endif
            msg_error_ = "No variable selected for add in the add step.";
            break;
          }

          action.second.erase(action.second.begin(),action.second.end());
          for(vector<int>::iterator it = newId.begin() ; it != newId.end(); it++)
            updateR(*it,signC,action);
        }
        else
          action=make_pair(false,dropId);

        addCmax(Cmax);

        //compute the inverse of G
        computeGi1(Gi1,signC);

        //compute Aa
        Aa = 1/sqrt(Gi1.sum());

        //compute w
        w = Gi1*signC*Aa;

        //compute equiangular vector
        u = Xi_*w;

        //computation of gamma hat
        //if the number of active variable is equal to the max number authorized, we don't search a new index
        if( nbActiveVariable_ == min(n_-1, p_-nbIgnoreVariable_) )
           gam = Cmax/Aa;
        else
        {
          //computation of a
          a = X_.transpose()*u;

          //computation of gamma hat
          gam = computeGamHat(Aa,a,Cmax);
        }

        //computation of gamma tilde
        gammaTilde = computeGamTilde(w,dropId);

//        cout<<"gamhat "<<gam<<"  gamtilde "<<gammaTilde<<std::endl;
//        for(int i=1;i<=path_.lastState().size();i++)
//          stk_cout<<path_.lastState()[i].first<<"       ";
//        stk_cout<<endl;
//        for(int i=1;i<=path_.lastState().size();i++)
//          stk_cout<<path_.lastState()[i].second<<" ";
//        stk_cout<<endl;
//        stk_cout<<signC;
//        stk_cout<<w;
//        stk_cout<<signC.size()<<"   "<<path_.lastState().size()<<"   "<<nbActiveVariable_<<std::endl;

        if( gammaTilde < gam )
        {
          gam = gammaTilde;
          isAddCase = false;
//          cout<<"to drop: ";
          for(int i = 0; i<(int) dropId.size(); i++)
          {
            nbActiveVariable_--;
//            cout<<activeVariables_[dropId[i]]<<"  ";
          }
//          cout<<endl;

//          cout<<action.second.size()<<" addcancel: ";
//          for(int i = 0; i < (int) action.second.size(); i++)
//          {
//            cout<<activeVariables_[action.second[i]]<<"  ";
//          }
//          cout<<endl;
        }
        else
          isAddCase = true;

        //update beta
        updateBeta(gam,w,action,isAddCase,dropId);

        //update of c
        c_ -= (X_.transpose()*u)*gam;

        //drop situation
        if(!isAddCase)
          dropStep(dropId,signC);

        //path_.states(step_).printCoeff();
      }
    }
    Real t1 = Chrono::elapsed();

#ifdef VERBOSE
      cout<<endl<<"Algorithm finished in "<<t1<<"s"<<endl;
      cout<<"Number of steps: "<<step_<<endl;
      cout<<"Number of active variables: "<<nbActiveVariable_<<endl;
#endif

  }


  /*
   * predict the path for a ratio fraction = l1norm/l1normmax
   * @param X new data for predict the response
   * @param fraction real between 0 and 1 .
   * @return predicted response
   */
  void Lars::predict(STK::CArrayXX const& X, STK::Real fraction, STK::CVectorX &yPred)
  {
    yPred = mu_;

    //fraction = 0 : all coefficients are equal to 0
    if(fraction == 0)
      return ;

    //fraction = 1 : coefficients of the last step
    if (fraction == 1)
    {
      int lastStep = path_.size()-1;//stocké dans un vector index à 0
      int nbVar = path_.lastState().sizeRows();


      for(int i = 1; i <= yPred.sizeRowsImpl(); i++)
        for(int j = 1; j <= nbVar; j++)
          yPred[i] += X(i, varIdx(lastStep,j)) * coefficient(lastStep,j);

      return ;
    }

    //fraction >0 and <1
    Array2DVector<STK::Real> l1norm(path_.l1norm());

    fraction *= l1norm.back();

    int index = 1;
    while(l1norm[index] < fraction)
      index++;

    //compute coefficient
    STK::Array2DVector< pair<int,Real> > coeff(std::max(path_.states(index-2).size(),path_.states(index-1).size()));
    //coeff.move(computeCoefficients(path_.states(index-2),path_.states(index-1),path_.evolution(index-2),fraction));
    computeCoefficients(path_.states(index-2),path_.states(index-1),path_.evolution(index-2),fraction,coeff);

    for( int i = 1; i <= yPred.sizeRowsImpl(); i++)
      for( int j = 1; j <= coeff.sizeRows(); j++)
        yPred[i] += X(i, coeff[j].first) * coeff[j].second ;

    return ;
  }


}//end namespace
