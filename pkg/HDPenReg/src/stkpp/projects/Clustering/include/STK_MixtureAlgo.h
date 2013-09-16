/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2012  Serge Iovleff

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
 * Project:  stkpp::Clustering
 * created on: 16 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_IMixtureModel.h
 *  @brief In this file we define the classes for mixture algorithms.
 **/


#ifndef STK_MIXTUREALGO_H
#define STK_MIXTUREALGO_H

#include "STK_Clust_Util.h"

#include "../../Sdk/include/STK_IRunner.h"
#include "../../STKernel/include/STK_Real.h"

namespace STK
{
// forward declaration
class IMixtureModelBase;

/** @ingroup Clustering
 * Interface base class for the algorithms.
 * All algorithms are runners applying on a model instance given by pointer
 * and have to implement the run method.
 **/
class IMixtureAlgo : public IRunnerBase
{
  protected:
    /** default constructor */
    inline IMixtureAlgo() : IRunnerBase(), p_model_(0), nbIterMax_(Clust::maxIterShortRun), epsilon_(Clust::epsilonShortRun) {}
    /** Constructor. Instantiate the algorithm with a mixture model.
     *  @param p_model a pointer on an instance of a mixture model */
    inline IMixtureAlgo( IMixtureModelBase* p_model, int nbIterMax, int epsilon)
                       : IRunnerBase(), p_model_(p_model), nbIterMax_(nbIterMax), epsilon_(epsilon)
    {}
    /** Copy constructor.
     *  @param mixtureAlgo the mixture to copy */
    inline IMixtureAlgo( IMixtureAlgo const& mixtureAlgo) : IRunnerBase(mixtureAlgo)
                       , p_model_(mixtureAlgo.p_model_), nbIterMax_(mixtureAlgo.nbIterMax_), epsilon_(mixtureAlgo.epsilon_)
    {}

  public:
    /** destructor */
    inline virtual ~IMixtureAlgo() {}
    /** set a new model */
    inline void setModel(IMixtureModelBase* p_model) { p_model_ = p_model; }
    /** set the maximal number of iterations */
    inline void setNbIterMax(int nbIterMax) { nbIterMax_ = nbIterMax; }
    /** set the maximal number of iterations */
    inline void setEpsilon(int epsilon) { epsilon_ = epsilon; }

  protected:
    /** pointer on the mixture model */
    IMixtureModelBase* p_model_;
    /** number of iterations of the algorithm */
    int nbIterMax_;
    /** tolerance of the algorithm. */
    Real epsilon_;
};

/** @ingroup Clustering
 *  Implementation of the EM algorithm.
 **/
class EMAlgo: public IMixtureAlgo
{
  public:
    /** default constructor */
    inline EMAlgo() : IMixtureAlgo() {}
    /**  Instantiate the algorithm with a mixture model.
     * @param p_model a pointer on an instance of a mixture model */
    inline EMAlgo( IMixtureModelBase* p_model, int nbIterMax, int epsilon)
                 : IMixtureAlgo(p_model, nbIterMax, epsilon) {}
    /** Copy constructor.
     *  @param emAlgo the mixture to copy */
    inline EMAlgo( EMAlgo const& emAlgo) : IMixtureAlgo(emAlgo) {}
    /** destructor */
    inline virtual ~EMAlgo(){}
    /** run the algorithm on the model calling eStep and mStep
     *  until the maximal number of iteration is reached or the variation
     *  of the lnLikelihood is less than epsilon.
     * @return @c true if no error occur, @c false otherwise*/
    virtual bool run();
};

/** @ingroup Clustering
 *  Implementation of the CEM algorithm.
 **/
class CEMAlgo: public IMixtureAlgo
{
  public:
    /** default constructor */
    inline CEMAlgo() : IMixtureAlgo() {}
    /**  Instantiate the algorithm with a mixture model.
     * @param p_model a pointer on an instance of a mixture model */
    inline CEMAlgo( IMixtureModelBase* p_model, int nbIterMax, int epsilon)
                  : IMixtureAlgo(p_model, nbIterMax, epsilon)
    {}
    /** destructor */
    inline virtual ~CEMAlgo(){}
    /** run the algorithm  on the model calling ceStep and mStep
     *  until the maximal number of iteration is reached or the variation
     *  of the lnLikelihood is less than epsilon..
     * @return @c true if no error occur, @c false otherwise*/
    virtual bool run();
};

/** @ingroup Clustering
 *  Implementation of the SEM algorithm.
 **/
class SEMAlgo: public IMixtureAlgo
{
  public:
    /** default constructor */
    inline SEMAlgo() : IMixtureAlgo() {}
    /** Instantiate the algorithm with a mixture model.
     *  @param p_model a pointer on an instance of a mixture model */
    inline SEMAlgo( IMixtureModelBase* p_model, int nbIterMax, int epsilon)
                  : IMixtureAlgo(p_model, nbIterMax, epsilon)
    {}
    /** destructor */
    inline virtual ~SEMAlgo(){}
    /** run the algorithm on the model calling SeStep
     *  until the maximal number of iteration is reached.
     *  @return @c true if no error occur, @c false otherwise. */
    virtual bool run();
};
} // namespace STK

#endif /* STK_MIXTUREALGO_H */
