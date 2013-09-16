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
 **/

/** @file STK_IMixtureModelBase.h
 *  @brief In this file we define the abstract base class for mixture models.
 **/

#ifndef STK_IMODELMIXTUREBASE_H
#define STK_IMODELMIXTUREBASE_H

#include "../../Arrays/include/STK_Array2D.h"
#include "../../StatModels/include/STK_IModelBase.h"

namespace STK
{

/** @brief base class for the  mixture model.
 * In this interface we assume there is an underline generative model that will
 * be estimated using either an EM, SEM or CEM algorithm.
 *
 * All mixture parameters: proportions, Tik, Zi and components are accessed by
 * pointer and can be set using the method
 * @code
 *   void setMixtureParameters(Array2D<Real>* p_prop, Array2D<Real>* p_tik, Array2DVector<int>* p_zi);
 * @endcode
 * to this class so that they can be used in a composed model.
 *
 * They also can be created using the method
 * @code
 *   void createMixtureParameters();
 * @endcode
 *
 * The pure virtual function to implement in derived class are
 * @code
 *   virtual IMixtureModelBase* create() const = 0;
 *   virtual IMixtureModelBase* clone() const = 0;
 *   virtual bool randomInit() =0;
 *   virtual void mStep() = 0;
 *   virtual Real componentProbability(int i, int k) = 0;
 * @endcode
 *
 * @note the range of the samples and the number of cluster have to be set
 * before any use of this class.
 * @note the proportions are computed in the mStep() pure virtual method in
 * order to allow equal (or user-fixed) proportions.
 */
class IMixtureModelBase : public IModelBase
{
  protected:
    /** default constructor */
    IMixtureModelBase();
    /** copy constructor. If the pointer on the mixture parameters are not zero
     *  then they are cloned.
     *  @param model the model to clone
     **/
    IMixtureModelBase( IMixtureModelBase const& model);

  public:
    /** destructor */
    virtual ~IMixtureModelBase();

    /** create pattern */
    virtual IMixtureModelBase* create() const = 0;
    /** clone pattern */
    virtual IMixtureModelBase* clone() const = 0;
    /** initialize randomly the parameters of the components of the model */
    virtual bool randomInit() =0;
    /** estimate the proportions and the parameters of the components of the
     *  model given the current tik/zi mixture parameters values.
     **/
    virtual void mStep() = 0;
    /** @return the value of the probability of the sample sample in  teh component k.
     *  @param index of the sample
     *  @param k index of the component
     **/
    virtual Real componentProbability(int i, int k) = 0;

    /** write the parameters of the model in the stream os. */
    virtual void writeParameters(std::ostream& os) const {};
    /** initialize randomly the labels zi of the model */
    void classInit();
    /** initialize randomly the posterior probabilities tik of the model */
    void fuzziInit();
    /** compute tik and zi, replace tik by hard classification. */
    void ceStep();
    /** compute tik and simulate zi, replace tik by hard classification.  */
    void seStep();
    /** compute tik, default implementation. */
    void eStep();
    /** Compute zi using the Map estimator, default implementation. */
    void mapStep();
    /** compute the ln-likelihood of the mixture model. */
    void computeLnLikelihood();

    /** @return the number of cluster. */
    inline int nbCluster() const { return nbCluster_;}
    /** @return the proportions of each mixtures */
    Array2DPoint<Real> const* p_prop() const { return p_prop_;}
    /** @return a constant pointer on the tik probabilities */
    Array2D<Real> const* p_tik() const { return p_tik_;}
    /** @return a constant pointer on the zi labels */
    Array2DVector<int> const* p_zi() const { return p_zi_;}

    /** set the range of the samples */
    inline void setRangeSamples( Range rangeSamples) { rangeSamples_ = rangeSamples;}
    /** set the proportions */
    inline void setNbCluster(int nbCluster) { nbCluster_ = nbCluster;}

    /** set the parameters of the  mixture model.
     *  @param p_prop pointer on the proportion of the mixture model
     *  @param p_tik pointer on the posterior probabilities
     *  @param p_zi pointer on the class labels
     * */
    void setMixtureParameters(Array2DPoint<Real>* p_prop, Array2D<Real>* p_tik, Array2DVector<int>* p_zi);

    /** Create the parameters of the  mixture model. */
    void createMixtureParameters();

  protected:
    /** range of the samples. */
    Range rangeSamples_;
    /** number of cluster. */
    int nbCluster_;

    /** The proportions of each mixtures */
    Array2DPoint<Real>* p_prop_;
    /** The tik probabilities */
    Array2D<Real>* p_tik_;
    /** The zik class label */
    Array2DVector<int>* p_zi_;

  private:
    /** Boolean checking if the mixture parameters have been created or set by the
     *  end-user*/
    bool isParametersCreated_;
    /** create the proportions and initialize them with equal values*/
    void createProp();
    /** create the tik probabilities array and initialize them with equal values*/
    void createTik();
    /** create the zi labels array and initialize them with equal values */
    void createZi();
    /** replace tik by zik, the indicator variable of the zi */
    void computeZik();
};

} // namespace SDTK

#endif /* IMODEL_H_ */
