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
 * created on: 6 juin 2013
 * Author:   Quentin Grimonprez
 **/

/** @file Cvlars.cpp
 *  @brief In this file, implementation of the methods of @c Cvlars .
 **/

#include "Cvlars.h"
#include <cstdlib>
#include <algorithm>

namespace HD
{
  /*
   * Constructor with no index ( it will be a sequence from 0 to 1 by 0.01)
   * @param X matrix of data, a row=a individual
   * @param y response
   * @param k number of folds
   * @param maxSteps number of maximum step to do
   * @param eps epsilon (for 0)
   */
  Cvlars::Cvlars(STK::CArrayXX const& X, STK::CVectorX const& y, int k, int maxSteps, bool intercept, STK::Real eps)
                : p_X_(&X)
                , p_y_(&y)
                , partition_(X.sizeRows())
                , sizePartition_(k,0)
                , index_(101)
                , lambdaMode_(false)
                , residuals_(101,k)
                , cv_(101,0)
                , cvError_(101)
                , k_(k)
                , n_(X.sizeRows())
                , p_(X.sizeCols())
                , maxSteps_(maxSteps)
                , eps_(eps)
                , intercept_(intercept)
  {
    //no index given, we created a sequence of inde between 0 and 1
    for(int i = 0; i<101; i++)
      index_[i] = (double) i/100;

    //create the partition
    partition();
  }

  /*
   * Constructor
   * @param X matrix of data, a row=a individual
   * @param y response
   * @param k number of folds
   * @param index vector with real between 0 and 1 (ratio (norm coefficient)/max(norm coefficient) for which we compute the prediction error)
   * @param maxSteps number of maximum step to do
   * @param eps epsilon (for 0)
   */
  Cvlars::Cvlars(STK::CArrayXX const& X, STK::CVectorX const& y, int k, std::vector<double> const& index, bool lambdaMode,int maxSteps, bool intercept, STK::Real eps)
                : p_X_(&X)
                , p_y_(&y)
                , partition_(X.sizeRows())
                , sizePartition_(k,0)
                , index_(index)
                , lambdaMode_(lambdaMode)
                , residuals_(index.size(),k)
                , cv_(index.size(),0)
                , cvError_(index.size())
                , k_(k)
                , n_(X.sizeRows())
                , p_(X.sizeCols())
                , maxSteps_(maxSteps)
                , eps_(eps)
                , intercept_(intercept)
  {
    //create the partition
    partition();
  }

  /*
   * create a random partition in k folds
   */
  void Cvlars::partition()
  {
    //fill the container with the index of folds
    for(int i = 0 ; i< n_ ;i++)
    {
      partition_[i] = i%k_;
      sizePartition_[i%k_]++;
    }


    /*int pos=1;
    for(int i = 0; i<k_; i++)
    {
      //foldRange_[i] = STK::Range(pos,sizePartition_);
      pos += sizePartition_[i];
    }*/

    //make a random rearrangement
    srand(time(NULL));
    random_shuffle(partition_.begin(),partition_.end());
  }

  void Cvlars::setPartition(std::vector<int> const& partition)
  {
    partition_ = partition;
    sizePartition_.resize(k_);
    for(int i = 0; i < k_; i++)
      sizePartition_[i] = 0;

    for(int i = 0; i < n_; i++)
      sizePartition_[partition_[i]]++;

  }

  /*
   * run the cross validation
   */
  void Cvlars::run()
  {
    //search the first and last fold with the same size
    std::vector<int> startIndex(1,0),endIndex(1,k_-1);
    int k = 0;
    for(int i = 1; i < k_; i++)
    {
      if(sizePartition_[i]!= sizePartition_[startIndex[k]])
      {
        startIndex.push_back(i);
        endIndex[k] = i-1;
        endIndex.push_back(k_-1);
        k++;
      }
    }

    //run for each size of fold
    for(int i = 0; i < (int) startIndex.size(); i++)
      subrun(startIndex[i],endIndex[i]);

    // compute mean prediction error for each index
    STK::CVectorX one(k_,1);
    cv_ = (residuals_ * one) / k_;

    // compute mean standard deviation of cv_ for each index
    for(int i = 1; i <= (int) index_.size(); i++)
      residuals_.row(i) -= cv_[i];

    residuals_ = residuals_.square();
    cvError_ = (residuals_ * one)/(k_-1)/k_;
    cvError_ = cvError_.sqrt();

  }

  /*
   * run cross validation for folds from idxStartFold to idxEndFold
   * @param idxStartFold index of the first fold
   * @param idxEndFold index of the last fold
   */
  void Cvlars::subrun(int idxStartFold,int idxEndFold)
  {
    //create test and control container
    STK::CArrayXX XControl( n_ - sizePartition_[idxStartFold], p_);
    STK::CVectorX yControl( n_ - sizePartition_[idxStartFold] );
    STK::CArrayXX XTest(sizePartition_[idxStartFold], p_);
    STK::CVectorX yTest(sizePartition_[idxStartFold] );
    STK::CVectorX yPred(sizePartition_[idxStartFold] );

    for(int i = idxStartFold ; i <= idxEndFold ; i++)
    {
      //fill the container
      int index = 1;
      int index2 = 1;
      for(int j = 1; j <= n_; j++)
      {
        if(partition_[j-1] != i)
        {
          yControl[index] = (*p_y_)[j];
          XControl.row(index)=p_X_->row(j);
          index++;
        }
        else
        {
          yTest[index2] = (*p_y_)[j];
          XTest.row(index2)=p_X_->row(j);
          index2++;
        }
      }

      //run lars on control data set
      HD::Lars lars(XControl,yControl,maxSteps_,intercept_,eps_);
      lars.run();

      for(int s = 1 ; s <= (int) index_.size(); s++)
      {
        //we compute the prediction of the y associated to XTest
        lars.predict(XTest,index_[s-1], lambdaMode_, yPred);

        //compute the residuals
        residuals_(s,i+1) = (yPred-yTest).square().sum()/sizePartition_[i];
      }
    }
  }

#ifdef _OPENMP
  void Cvlars::run2()
   {
     //search the first and last fold with the same size
     std::vector<int> startIndex(1,0),endIndex(1,k_-1);
     int k = 0;
     for(int i = 1; i < k_; i++)
     {
       if(sizePartition_[i]!= sizePartition_[startIndex[k]])
       {
         startIndex.push_back(i);
         endIndex[k] = i-1;
         endIndex.push_back(k_-1);
         k++;
       }
     }

     //run for each size of fold
     //create test and control container
     #pragma omp parallel
     {
       #pragma omp for schedule(dynamic,1)
       for(int i = 0; i < k_ ; i++)
       {

         STK::CArrayXX XControl( n_ - sizePartition_[i], p_);
         STK::CVectorX yControl( n_ - sizePartition_[i] );
         STK::CArrayXX XTest(sizePartition_[i], p_);
         STK::CVectorX yTest(sizePartition_[i] );
         STK::CVectorX yPred(sizePartition_[i] );

         //fill the container
         int index = 1;
         int index2 = 1;
         for(int j = 1; j <= n_; j++)
         {
           if(partition_[j-1] != i)
           {
             yControl[index] = (*p_y_)[j];
             XControl.row(index)=p_X_->row(j);
             index++;
           }
           else
           {
             yTest[index2] = (*p_y_)[j];
             XTest.row(index2)=p_X_->row(j);
             index2++;
           }
         }

         //run lars on control data set
         HD::Lars lars(XControl,yControl,maxSteps_,intercept_,eps_);
         lars.run();

         for(int s = 1 ; s <= (int) index_.size(); s++)
         {
           //we compute the prediction of the y associated to XTest
           lars.predict(XTest,index_[s-1], lambdaMode_, yPred);

           //compute the residuals
           residuals_(s,i+1) = (yPred-yTest).square().sum()/sizePartition_[i];
         }
       }
     }//end parallel

     // compute mean prediction error for each index
     STK::CVectorX one(k_,1);
     cv_ = (residuals_ * one) / k_;

     // compute mean standard deviation of cv_ for each index
     for(int i = 1; i <= (int) index_.size(); i++)
       residuals_.row(i) -= cv_[i];
     residuals_ = residuals_.square();
     cvError_ = (residuals_ * one)/(k_-1)/k_;
     cvError_ = cvError_.sqrt();

   }
#else
  void Cvlars::run2()
   {
     run();
   }
#endif

}//end namespace HD
