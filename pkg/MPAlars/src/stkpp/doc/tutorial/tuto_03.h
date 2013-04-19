/* Contents of the tutorials
 **/
namespace STK
{

/** @page TutorialStatModel How to define a statistical model
 *
 *  @section StatModelIntroduction The IMultiStatModel interface
 *  All the multivariate statistical models should inherit from the
 *  templated interface class
 *  @code
 *  template <class Array, class MultivariateLaw> class IMultiStatModel;
 *  @endcode
 *  The first template parameter Array represent the type of a two dimensional
 *  Array deriving from ITContainer. It should be an Array2D or a
 *  CArray container. The second template parameter is the type
 *  of a Multivariate law deriving from the interface class Law::IMultiLaw or,
 *  in case of a joint probability, from the Interface class Law::JointProbability.
 *
 *  The class IMultiStatModel inherit from:
 *  @li IModelBase, which furnish the following accessors
 *  @code
 *    int const& nbSample() const;
 *    Real lnNbSample() const;
 *    int const& nbVar() const;
 *    Real lnLikelihood() const;
 *    Real likelihood() const;
 *    int const& nbFreeParameter();
 *  @endcode
 *  The number of samples and variables are determined by the size of the Array
 *  while the other fields have to be computed and updated by derived classes.
 *
 *  @li IRunnerConst<typename Array::Col>, which ask for implementing the
 *  following virtual functions
 *  @code
 *    bool run() =0;
 *    bool run(Array::Col const& weights) =0;
 *    void update() {};
 *  @endcode
 *
 *  The first @c run() method is used when the user of the class want to estimate
 *  the parameters of the model by, for example, maximizing the log-likelihhod.
 *  The second @c run(weights) methods is used when the user of the class want
 *  to estimate the parameters of the model by maximizing the weighted log-likelihood.
 *  @c udpate() is to be overloaded if when assigning a data set there is some
 *  task to process.
 *
 *  @section StatModelArray Prerequisites for the Array
 *  The Array have to provide two typedef @c Row and @c Col and two
 *  accessors @c row and @c col which can be used in the following form
 *  @code
 *    Array data;
 *    data.col(j); // return a typename Array::Col column vector
 *    data.row(i); // return a typename Array::Row row vector
 *  @endcode
 *
 *  @section StatModelLaw Prerequisites for the probability law
 *
 *  @subsection StatModelMultiLaw The Law::IMultilaw class
 *
 *  The class Law::IMultilaw is a templated interface class
 *  @code
 *  template<class RowVector> class IMultiLaw;
 *  @endcode
 *  The typedef @c RowVector must be the same than the typedef @c Array::Row.
 *
 *  The class IMultLaw derive from the interface class ILawBase and have to
 *  implement the following pure virtual methods
 *  @code
 *    virtual Real pdf( RowVector const& x) const =0;
 *    virtual Real lpdf( RowVector const& x) const =0;
 *    virtual void rand( RowVector& x) const =0;
 *  @endcode
 *  It is a very general Interface that can be used when data are correlated.
 *
 *  @subsection StatModelJointLaw The Law::JointProbability class
 *  The Law::JointProbability is derived from IMultiLaw and is defined as
 *  @code
 *    template<class RowVector, class Law> class JointProbability;
 *  @endcode
 *  the template class RowVector is described in @ref StatModelMultiLaw.
 *  The template class Law, is expecting a class derived from the interface
 *  class Law::IUnivLaw.
 *
 *  For example the JointBernoulli class is defined as
 *  @code
 *  template<class RowVector>
 *   class JointBernoulli: public JointProbability<RowVector, Bernoulli>
 *  @endcode
 *  and allow to modelize a Binary data set when assuming that
 *  each column is a sample of an independent Bernoulli law.
 *  @sa Law::JointBernoulli
 *
 *  @section StatModelExample Example
 *  The JointBernoulliModel is a templated class defined as
 *  @code
 *    template <class Array> class JointBernoulliModel;
 *  @endcode
 *  and inherit from
 *  @code
 *   IMultiStatModel<Array, Law::JointBernoulli< typename Array::Row> >
 *  @endcode
 *  @sa JointBernoulliModel
 **/

} // namespace STK
