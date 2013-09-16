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
 * Project:  stkpp::Clustering
 * created on: 24 août 2013
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_MixtureInit.h
 *  @brief In this file we define the initialization methods.
 **/


#ifndef STK_MIXTUREINIT_H
#define STK_MIXTUREINIT_H


#include "../../Sdk/include/STK_IRunner.h"

namespace STK
{
// forward declaration
class IMixtureModelBase;

/** @ingroup Clustering
 *  Interface base class for the initializations.
 *  All derived class will apply on a model instance and have to implement the
 *  run method.
 **/
class IMixtureInit : public IRunnerBase
{
  protected:
    /** default constructor */
    inline IMixtureInit() : IRunnerBase(), p_model_(0) {}
    /** copy constructor.
     * @param init the initializing method to copy
     **/
    inline IMixtureInit(IMixtureInit const& init) : IRunnerBase(init), p_model_(init.p_model_) {}
    /**  Instantiate the algorithm with a mixture model.
     * @param p_model a pointer on an instance of a mixture model */
    inline IMixtureInit(IMixtureModelBase* p_model) : p_model_(p_model) {}

  public:
    /** destructor */
    inline virtual ~IMixtureInit() {}
    /** set a new model */
    inline void setModel(IMixtureModelBase* p_model) { p_model_ = p_model; }

  protected:
    /** pointer on the mixture model */
    IMixtureModelBase* p_model_;
};

/** @ingroup Clustering
 *  Implementation of the random initialization. This class will initialize the
 *  parameter by calling the randomInit() method of the model. */
class RandomInit: public IMixtureInit
{
  public:
    /** default constructor */
    inline RandomInit() : IMixtureInit() {}
    /**  Instantiate the algorithm with a mixture model.
     * @param p_model a pointer on an instance of a mixture model */
    inline RandomInit(IMixtureModelBase* p_model) : IMixtureInit(p_model) {}
    /** destructor */
    inline virtual ~RandomInit(){}
    /** run the initialization by calling the randomInit method of the model.
     * @return @c true if no error occur, @c false otherwise*/
    virtual bool run();
};

/** @ingroup Clustering
 *  Initialization by simulating a realization of the class labels zi accordingly
 *  to the initial proportions.
 **/
class ClassInit: public IMixtureInit
{
  public:
    /** default constructor */
    inline ClassInit() : IMixtureInit() {}
    /**  Instantiate the algorithm with a mixture model.
     * @param p_model a pointer on an instance of a mixture model */
    inline ClassInit(IMixtureModelBase* p_model) : IMixtureInit(p_model) {}
    /** destructor */
    inline virtual ~ClassInit(){}
    /** run the initialization by calling the classInit method of the model.
     * @return @c true if no error occur, @c false otherwise*/
    virtual bool run();
};

/** @ingroup Clustering
 *  Initialization by simulating the tik accordingly to the initial
 *  proportions. */
class FuzziInit: public IMixtureInit
{
  public:
    /** default constructor */
    inline FuzziInit() : IMixtureInit() {}
    /**  Instantiate the algorithm with a mixture model.
     * @param p_model a pointer on an instance of a mixture model */
    inline FuzziInit(IMixtureModelBase* p_model) : IMixtureInit(p_model) {}
    /** destructor */
    inline virtual ~FuzziInit(){}
    /** run the algorithm on the model calling E step and M step.
     * @return @c true if no error occur, @c false otherwise*/
    virtual bool run();
};

} // namespace STK

#endif /* STK_MIXTUREINIT_H */
