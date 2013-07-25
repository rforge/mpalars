#include "stkpp/projects/STKernel/include/STK_Real.h"
#include "emR.h"
#include "larsR.h"
#include "lassoModels/PenalizedModels.h"
#include "lassoModels/EM.h"
#include "lassoModels/LassoPenalty.h"
#include "lassoModels/LassoSolver.h"

#include <iostream>

using namespace Rcpp;
using namespace STK;
using namespace std;
using namespace HD;


struct LeastSquareMultiplicator
{
	STK::CVectorX operator()(STK::CVectorX const& x) const
	{
		STK::CVectorX a;
		a = ( p_data_->transpose() * ( (*p_data_) * x ) );

		return   a ;
	}

	LeastSquareMultiplicator(STK::CArrayXX const* p_data) : p_data_(p_data) {}

	STK::CArrayXX const* p_data_;
};


void copySTKVectorInSTDVector(STK::CVectorX const& stkvector, vector<STK::Real> &stdvector)
{
	stdvector.resize(stkvector.sizeRows());
	for(int i = 1; i <= stkvector.sizeRows(); i++)
		stdvector[i-1]=stkvector[i];	
}

void copySTKArray2DVectorInSTDVector(STK::Array2DVector<int> const& stkvector, vector<int> &stdvector)
{
	stdvector.resize(stkvector.sizeRows());
	for(int i = 1; i <= stkvector.sizeRows(); i++)
		stdvector[i-1]=stkvector[i];	
}

RcppExport SEXP EMlasso(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP lambda, SEXP maxStep, SEXP eps)
{
  double t1,t2;
  t1=clock();

  //convert parameters
  int p(as<int>(nbVar)), n(as<int>(nbIndiv)), maxStepC(as<int>(maxStep));
  Real epsC(as<Real>(eps));

  CArrayXX x(n,p);
  convertToArray(data,x);
  CVectorX y(n);
  convertToVector(response,nbIndiv,y);

  vector<STK::Real> lambdaC = as<vector<STK::Real> >(lambda);


  //center the data for intercept
  STK::CVectorX muX(p);
  STK::Real mu;

  mu = y.sum()/n;
  y -= mu;

  for(int j = 1; j <= p; j++)
  {
    muX[j] = x.col(j).sum()/n;
    for(int i=1; i <= n; i++)
      x(i,j) -= muX[j];
  }

  //Xty for LassoSolver
  STK::CVectorX Xty(p);
  Xty=x.transpose() * y;

  //if lambdaC[0]=-1, we have to generate the lambda sequence with the same way as the glmnet package
  if(lambdaC[0] == -1)
  {
    lambdaC.resize(100);
    int pos;
    //corrMax is the first lambda value in the lars sequence
    STK::Real corrMax =Xty.abs().maxElt(pos);
    STK::Real lambdaMinRatio, minLambda;
    if(n<p)
      lambdaMinRatio=0.01;
    else
      lambdaMinRatio=0.0001;

    minLambda = corrMax * lambdaMinRatio;
    STK::Real gapLambda = (log(corrMax)-log(minLambda))/99;
    lambdaC[99] = corrMax;
    for(int i = 99; i>0 ; i--)
      lambdaC[i-1] = exp(log(lambdaC[i]) - gapLambda); 	
  }
	


  //container for all the solution
  vector<vector<int> > pathSolutionIndex(lambdaC.size()+1);//no 0 in lambda, we had the LS solution in first position
  vector<vector<STK::Real> > pathSolutionCoefficient(lambdaC.size()+1);//no 0 in lambda, we had the LS solution in first position

  vector<int> indexTemp;
  vector<STK::Real> coefficientTemp;
  t2=clock();


  //find the Least square solution
  LeastSquareMultiplicator multLS(&x);

  CG<LeastSquareMultiplicator,STK::CVectorX> LSsolver(multLS,Xty,0,epsC);
  LSsolver.run();

  //create PenalizedModels
  PenalizedModels lasso(&x,y,LSsolver.x(),0,0,epsC);
  PenalizedModels* p_lasso(&lasso);

  //add the LS solution 
  copySTKVectorInSTDVector(lasso.currentBeta(), coefficientTemp);
  pathSolutionCoefficient[0]=coefficientTemp;
  copySTKArray2DVectorInSTDVector(lasso.currentSet(), indexTemp);
  pathSolutionIndex[0]=indexTemp;

  //create penalty
  LassoPenalty penalty(lambdaC[0],n,p);//don't run for lambda = 0
  LassoPenalty* const p_penalty(&penalty);

  //add the penalty
  lasso.setPenalty(p_penalty);

  //create CG
  //creation functor for CG
  LassoMultiplicator mult(lasso.p_currentData(),penalty.p_invPenalty(), penalty.p_sigma2());
  //create solver
  InitFunctor init(lasso.p_currentBeta());//pointer on the currentBeta_ of penalizedModels
  InitFunctor* p_init(&init);
  CG<LassoMultiplicator,STK::CVectorX,InitFunctor> gcsolver(mult,Xty,p_init,epsC);
  CG<LassoMultiplicator,STK::CVectorX,InitFunctor>* p_gcsolver(&gcsolver);

  //create solver for lasso
  LassoSolver lassosolver(lasso.p_currentData(),lasso.p_currentSet(),Xty,p_gcsolver,p_penalty);
  LassoSolver* p_lassosolver(&lassosolver);
  lasso.setSolver(p_lassosolver);

  //create EM
  EM algo(p,epsC);


  //run for all lambda
  for(int i=1; i<=lambdaC.size(); i++)
  {
     //change the lambda
     penalty.setLambda(lambdaC[i-1]);

     //run algorithm
     Chrono::start();
     algo.run(p_lasso);
     t1=Chrono::elapsed();

     //stock the new values
     copySTKVectorInSTDVector(lasso.currentBeta(), coefficientTemp);
     pathSolutionCoefficient[i]=coefficientTemp;
     copySTKArray2DVectorInSTDVector(lasso.currentSet(), indexTemp);
     pathSolutionIndex[i]=indexTemp;
  }


  t2=clock();


  return List::create(Named("index")=wrap(pathSolutionIndex),Named("coefficient")=wrap(pathSolutionCoefficient),Named("lambda")=wrap(lambdaC));

}

