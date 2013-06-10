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
 * Project:  stkpp::Model
 * created on: 22 juil. 2011
 * Purpose: define the class IMultiStatModel.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_IMultiStatModel.h
 *  @brief In this file we define the class IMultiStatModel.
 **/

#ifndef STK_IMULTISTATMODEL_H
#define STK_IMULTISTATMODEL_H

#include <cmath>

#include "STK_IModelBase.h"
#include "../../STKernel/include/STK_Macros.h"
#include "../../STatistiK/include/STK_Law_IMultiLaw.h"

namespace STK
{

/** @ingroup StatModels
 *  @brief Interface base class for the parameters of a multivariate model.
  */
class IMultiParameters
{
  protected:
    /** default constructor.*/
    IMultiParameters() {}

  public:
    /** Destructor */
    virtual ~IMultiParameters() {}
    /** clone the parameters
     *  @return a clone of the current parameters
     **/
    virtual IMultiParameters* clone() const = 0;
    /** resize the parameters.
     *  @param size the size of the parameters (range of the variables)
     **/
     virtual void resize( Range const& size) = 0;
};

/** @ingroup StatModels
 *  @brief Interface base class for all Multivariate Statistical Models.
 *
 *  A Statistical model, \f$ \mathcal{P}\f$, is a collection of multivariate
 *  probability distribution functions or probability density functions
 *  (collectively referred to as ''distributions'' for brevity).
 *  A parametric model is a collection of distributions, each of which is
 *  indexed by a unique finite-dimensional parameter:
 *  \f$\mathcal{P}=\{\mathbb{P}_{\theta} : \theta \in \Theta\}\f$, where
 *  \f$\theta\f$ is a parameter and \f$\Theta \subseteq \mathbb{R}^d\f$ is
 *  the feasible region of parameters, which is a subset of d-dimensional
 *  Euclidean space.
 *
 *  A statistical model may be used to describe the set of
 *  distributions from which one assumes that a particular data set is sampled.
 *  For example, if one assumes that data arise from a multivariate Gaussian
 *  distribution, then one has assumed a Gaussian model:
 *  \f$
 *    \mathcal{P}=\{\mathbb{P}(x; \mu, \Sigma) = \frac{1}{\sqrt{2 \pi |\Sigma|} }
 *    \exp\left\{ -\frac{1}{2}(x-\mu)'\Sigma^{-1}(x-\mu)\right\} : \mu \in \mathbb{R}^p, \Sigma > 0\}
 *  \f$.
 *
 *  From a computational point of view a statistical model is defined with
 *  the help of two elements
 *  - A data set where the number of samples is the number of rows and the number
 *  of variable is the number of columns. This data set is stored in a Container
 *  of type @c Array.
 *  - A set of parameters stored in a class of type @c Parameters.
 *
 *  Derived implementations of this interface have to implement the following
 *  pure methods:
 *  @code
 *    int computeNbFreeParameters() const;
 *    Real computeLnLikelihood( RowVector const& rowData) const = 0;
 *    virtual void computeParameters() = 0;
 *    virtual void computeParameters(ColVector const& weights) = 0;
 *  @endcode
 *
 *  @tparam Array can be any kind of vector for the data set. it should at
 *  least derive from the specialization of ITContainer for Arrays::vector_,
 *  @sa ITContainer.
 *  @tparam Parameters any structure encapsulating the parameters of the model.
 *
 *  @note this class can be a runner, in this case the parameters are estimated
 *  using either @c run() and @c run(weights) methods. This class can also be used
 *  as a "kitchen" providing tools, in particular if there is latent variables.
 *  @sa ILatentModel, MixtureModel.
 **/
template <class Array, class WColVector, class Parameters>
class IMultiStatModel : public IModelBase, public IRunnerUnsupervised<Array, WColVector>
{
  public:
    /** Type of the data contained in the container */
    typedef typename Array::Type Type;
    /** Type of the row vector of the container */
    typedef typename Array::Row RowVector;
    /** Type of the column vector of the container */
    typedef typename Array::Col ColVector;
    /** Type of the runner */
    typedef IRunnerUnsupervised<Array, WColVector> Runner;
    using Runner::p_data;

  protected:
    /** default constructor. */
    IMultiStatModel() : IModelBase(), Runner(), p_param_(0) {}
    /** Constructor with data set. */
    IMultiStatModel(Array const& data) : IModelBase(), Runner(data), p_param_(0)
    { initializeModel();}

    /** Constructor with a ptr on the data set. */
    IMultiStatModel(Array const* p_data) : IModelBase(), Runner(p_data), p_param_(0)
    { initializeModel();}
    /** Copy constructor.
     *  @param model the model to copy
     **/
    IMultiStatModel( IMultiStatModel const& model)
                   : IModelBase(model), Runner(model), p_param_(model.p_param_->clone())
    {}

  public:
    /** destructor */
    virtual ~IMultiStatModel() { if (p_param_) delete p_param_;}
    /** @return the pointer on the parameters */
    inline Parameters* const p_param() const { return p_param_;}
    /** @param p_param the pointer on the parameters to set */
    inline void setParameters(Parameters* p_param) { p_param_ = p_param;}
    /** @return the pointer on the parameters and release them from teh class. */
    inline Parameters* release() { Parameters* p =p_param_; p_param_=0; return p;}
    /** compute the log Likelihood of the statistical model. */
    Real computeLnLikelihood() const
    {
      Real sum = 0.0;
      for (int i= p_data()->firstIdxRows(); i<= p_data()->lastIdxRows(); i++)
      { sum += computeLnLikelihood(p_data()->row(i));}
      return(sum);
    }
    /** Estimate the parameters of the model and update the
     **/
    virtual bool run()
    {
      if (!p_data())
      { this->msg_error_ = STKERROR_NO_ARG(IMultiStatModel::run,data have not be set);
        return false;
      }
      if (!p_param())
      { this->msg_error_ = STKERROR_NO_ARG(IMultiStatModel::run(weights),parameters have not be set);
        return false;
      }
      try
      {
        // compute parameters
        computeParameters();
        // compute log-likelihood
        this->setLnLikelihood(computeLnLikelihood());
        // set the number of free parameters
        this->setNbFreeParameters(computeNbFreeParameters());
      }
      catch (Exception const& e)
      { this->msg_error_ = e.error(); return false;}
      return true;
    }
    /** compute the weighted empirical probability of success based on the observed
     *  variables. The NA values are discarded.
     *  @param weights the weights of the observations
     **/
    virtual bool run(WColVector const& weights)
    {
      if (!p_data())
      { this->msg_error_ = STKERROR_NO_ARG(IMultiStatModel::run(weights),data have not be set);
        return false;
      }
      if (!p_param())
      { this->msg_error_ = STKERROR_NO_ARG(IMultiStatModel::run(weights),parameters have not be set);
        return false;
      }
      try
      {
        // compute weighted parameters
        computeParameters(weights);
        // compute log-likelihood
        this->setLnLikelihood(computeLnLikelihood());
        // set the number of free parameters
        this->setNbFreeParameters(computeNbFreeParameters());
      }
      catch (Exception const& e)
      { this->msg_error_ = e.error(); return false;}
      return true;
    }
    /** compute the number of free parameters */
    virtual int computeNbFreeParameters() const =0;
    /** compute the log Likelihood of an observation. */
    virtual Real computeLnLikelihood( RowVector const& rowData) const = 0;

  protected:
    /** Pointer on the parameters of the model. */
    Parameters* p_param_;
    /** This method is called if the user set a new data set
     *  @sa IRunnerUnsupervised::setData */
    virtual void update()
    {
    }
    /** compute the parameters */
    virtual void computeParameters() = 0;
    /** compute the weighted parameters */
    virtual void computeParameters(WColVector const& weights) = 0;

  private:
    /** initialize the model and the parameters */
    void initializeModel()
    {
      if (p_data())
      { this->initialize(p_data()->sizeRows(), p_data()->sizeCols());
        if (!p_param())
        { p_param_ = new Parameters(p_data()->cols());}
        else p_param_->resize(p_data()->cols());
      }
    }
};

} // namespace STK

#endif /* STK_IMULTISTATMODEL_H */
