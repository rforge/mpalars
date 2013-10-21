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
 *  @brief In this file we define the main class for mixture models.
 **/


#ifndef STK_IMIXTUREMODEL_H
#define STK_IMIXTUREMODEL_H

#include "STK_IMixtureModelBase.h"
#include "STK_MixtureTraits.h"
#include "../../Arrays/include/STK_Array1D.h"

namespace STK
{
/**@ingroup Clustering
 * Main interface class for mixture models.
 * At this level we add the array of Components and a
 * pointer on the data set. The components are created in this class.
 *
 * @tparam Array can be any kind of container for the observable variables.
 * It should at least derive from ITContainer and provide an
 * access to a single row. @sa ITContainer, ICArray, IArray2DBase
 *
 * @tparam Components any structure encapsulating the components
 * and deriving from IMixtureComponent.
 * @sa IMixtureComponent, IMultiStatModel, IMultiParameters
 **/
template <class Derived >
class IMixtureModel : public IRecursiveTemplate<Derived>, public IMixtureModelBase
{
  public:
    typedef typename hidden::MixtureTraits<Derived>::Array Array;
    typedef typename hidden::MixtureTraits<Derived>::Component Component;
    typedef typename hidden::MixtureTraits<Derived>::Parameters Parameters;
    using IMixtureModelBase::p_tik;

  protected:
    /** Default constructor */
    IMixtureModel( int nbCluster)
                 : IMixtureModelBase(nbCluster), p_data_(0), components_(nbCluster, 0)
    {
      for (int k= components_.firstIdx(); k <= components_.lastIdx(); ++k)
      { components_[k] = new Component(p_data_);}
    }
    /** Constructor with data set.
     *  @param nbCluster the number of cluster
     *  @param data a reference on the data set
     **/
    IMixtureModel( int nbCluster, Array const& data)
                 : IMixtureModelBase(nbCluster), p_data_(&data), components_(nbCluster, 0)
    {
      this->initialize(p_data_->rows().size(), p_data_->cols().size());
      for (int k= components_.firstIdx(); k <= components_.lastIdx(); ++k)
      { components_[k] = new Component(p_data_);}
    }
    /** Constructor with pointer on the data set.
     *  @param nbCluster the numebr of cluster
     *  @param p_data a pointer on the data set
     **/
    IMixtureModel( int nbCluster, Array const* p_data)
                 : IMixtureModelBase(nbCluster), p_data_(p_data), components_(nbCluster, 0)
    {
      if (p_data_)
      {
        this->initialize(p_data_->rows().size(), p_data_->cols().size());
      }
      for (int k= components_.firstIdx(); k <= components_.lastIdx(); ++k)
      { components_[k] = new Component(p_data_);}
    }
    /** copy constructor.
     *  Call the clone method of the Components class.
     *  @param model the model to copy
     **/
    IMixtureModel( IMixtureModel const& model)
                 : IMixtureModelBase(model), p_data_(model.p_data_), components_(model.components_)
    {
      for (int k= components_.firstIdx(); k <= components_.lastIdx(); ++k)
      { components_[k] = model.components_[k]->clone();}
    }
  public:
    /** destructor */
    virtual ~IMixtureModel()
    {
      for (int k= components_.firstIdx(); k <= components_.lastIdx(); ++k)
      { if (components_[k]) delete components_[k];}
    }
    /** clone pattern */
    IMixtureModel* clone() const
    { return IRecursiveTemplate<Derived>::clone();}
    /** create pattern */
    IMixtureModel* create() const
    {
      IMixtureModel* p_model =new Derived(this->nbCluster(), this->p_data_);
      if (this->isParametersCreated())
      { p_model->createMixtureParameters();}
      else
      { p_model->setMixtureParameters(this->p_prop_, this->p_tik_, this->p_zi_);}
      p_model->initializeModel();
      return p_model;
    }
    /** @return a pointer on the current data set */
    Array const* p_data() const { return p_data_;};
    /** @return the array with the components */
    Array1D< Component* > const& components() const { return components_;}
    /** @return a constant reference on the k-th component */
    Component* const& components(int k) const { return components_[k];}
    /** @brief Initialize the model before its first use.
     *  This function can be overloaded in derived class for initialization of
     *  the specific model parameters. It should be called prior to any used of
     *  the class. In this interface, the @c initializeModel() method
     *  - set the number of variables of the mixture model
     *  - set the range of the samples in the base class
     *  - call the base class IMixtureModelBase::initializeModel() method.
     **/
    virtual void initializeModel()
    {
      if (!p_data_)
      { STKRUNTIME_ERROR_NO_ARG(IMixtureModel,data set is not set);}
      this->initialize(p_data_->rows().size(), p_data_->cols().size());
      IMixtureModelBase::initializeModel();
    }
    /** set a new data set
     *  @param data the data set to set*/
    virtual void setData(Array const& data)
    {
      p_data_ = &data;
      this->initialize(p_data_->rows().size(), p_data_->cols().size());
      // update components
      for (int k= components_.firstIdx(); k <= components_.lastIdx(); ++k)
      { components_[k]->setData(p_data_);}
    }
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i index of the sample
     *  @param k index of the component
     **/
    virtual Real lnComponentProbability(int i, int k)
    { return components_[k]->computeLnLikelihood(p_data_->row(i));}

    /** Call static method mStep() implemented by end-user */
    virtual void mStep()
    { computeProportions();
      MixtureModelImpl< Array, Parameters, Component >::mStep(components_, p_tik());
    }
    /** use the default static method randomInit() implemented by end-users
     *  get an initial value of the parameters and compute tik and zi.
     **/
    virtual void randomInit()
    { MixtureModelImpl< Array, Parameters, Component >::randomInit(components_);
      eStep();
    }
    /** use the default static method initializeStep() for a first initialization
     *  of the parameters using tik values.
     **/
    virtual void initializeStep()
    { MixtureModelImpl< Array, Parameters, Component >::initializeStep(components_,  p_tik());}

  protected:
    /** pointer on the data set */
    Array const* p_data_;
    /** Array of the components of the mixture model */
    Array1D< Component* > components_;
};

/**@ingroup Clustering
 * Utility interface class for mixture models with fixed proportions.
 * In this class we overload the computeProportions() method
 * defined in the Interface base class IMixtureModelBase and let them
 * unchanged.
 **/
template <class Derived>
class IMixtureModelFixedProp : public IMixtureModel<Derived >
{
  protected:
    typedef typename hidden::MixtureTraits<Derived>::Array Array;
    typedef IMixtureModel<Derived> Base;
    /** Default constructor */
    IMixtureModelFixedProp( int nbCluster): Base(nbCluster)
    {}
    /** Constructor with data set.
     *  @param nbCluster the number of cluster
     *  @param data a reference on the data set
     **/
    IMixtureModelFixedProp( int nbCluster, Array const& data)
                          : Base(nbCluster, data)
    {}
    /** Constructor with pointer on the data set.
     *  @param nbCluster the numebr of cluster
     *  @param p_data a pointer on the data set
     **/
    IMixtureModelFixedProp( int nbCluster, Array const* p_data)
                          : Base(nbCluster, p_data)
    { }
    /** copy constructor.
     *  Call the clone method of the Components class.
     *  @param model the model to copy
     **/
    IMixtureModelFixedProp( IMixtureModelFixedProp const& model)
                 : Base(model)
    {}
    /** destructor */
    ~IMixtureModelFixedProp() {}
    /** overloading of the computePropotions() method.
     * Let them initialized to 1/K. */
    virtual void computeProportions() {}

};

} // namespace STK

#endif /* STK_IMIXTUREMODEL_H */
