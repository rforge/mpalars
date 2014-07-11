#include "stkpp/projects/STKernel/include/STK_Real.h"
#include "emR.h"
#include "larsR.h"
#include "HDPenReg/lassoModels/Lasso.h"
#include "HDPenReg/lassoModels/FusedLasso.h"
#include "HDPenReg/lassoModels/LogisticFusedLasso.h"
#include "HDPenReg/lassoModels/LogisticLasso.h"
#include "HDPenReg/lassoModels/EM.h"
#include "HDPenReg/lassoModels/CV.h"
#include "HDPenReg/lassoModels/CVLasso.h"
#include "HDPenReg/lassoModels/CVFusedLasso.h"

#include <iostream>

using namespace Rcpp;
using namespace STK;
using namespace std;
using namespace HD;



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

RcppExport SEXP EMlasso(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP lambda, SEXP intercept, SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG)
{
  //double t1,t2;
  //t1=clock();

  //convert parameters
  int p(as<int>(nbVar)), n(as<int>(nbIndiv)), maxStepC(as<int>(maxStep)), burnC(as<int>(burn));
  Real epsC(as<Real>(eps)), thresholdC(as<Real>(threshold)), epsCGC(as<Real>(epsCG));
  bool interceptC=as<bool>(intercept);

  CArrayXX x(n,p);
  convertToArray(data,x);
  CVectorX y(n);
  convertToVector(response,nbIndiv,y);

  vector<STK::Real> lambdaC = as<vector<STK::Real> >(lambda);

  //center the data for intercept
  STK::CVectorX muX(p);
  STK::Real mu = 0;

  if(interceptC)
  {
    mu = y.sum()/n;
    y -= mu;

    for(int j = 1; j <= p; j++)
    {
      muX[j] = x.col(j).sum()/n;
      for(int i=1; i <= n; i++)
        x(i,j) -= muX[j];
    }
  }


  //if lambdaC[0]=-1, we have to generate the lambda sequence with the same way as the glmnet package
  if(lambdaC[0] == -1)
  {
    STK::CVectorX Xty(p);
    Xty=x.transpose() * y;

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
  vector<vector<int> > pathSolutionIndex(lambdaC.size());
  vector<vector<STK::Real> > pathSolutionCoefficient(lambdaC.size());
  vector<int> step(lambdaC.size());
  vector<int> indexTemp;
  vector<STK::Real> coefficientTemp;

  //t2=clock();
  //create EM
  EM algo(maxStepC,burnC,epsC);
  Lasso lasso( &x, &y, lambdaC[0], thresholdC, epsCGC);

  //run for all lambda
  STK::CVectorX betaTemp(p);
  for(int i = 0; i < (int) lambdaC.size(); i++)
  {
    //change the lambda
    lasso.setLambda(lambdaC[i]);

    //run algorithm
    //Chrono::start();
    algo.run(&lasso);
    //t1=Chrono::elapsed();

    //stock the new values
    copySTKVectorInSTDVector(lasso.currentBeta(), coefficientTemp);
    pathSolutionCoefficient[i]=coefficientTemp;
    copySTKArray2DVectorInSTDVector(lasso.currentSet(), indexTemp);
    pathSolutionIndex[i]=indexTemp;
    step[i]=algo.step();

    //if there is 0 non-zeros coefficients, we stop at this values of lambda
    if( (coefficientTemp.size()==0))
    {
      lambdaC.erase(lambdaC.begin()+i+1,lambdaC.end());
      step.erase(step.begin()+i+1,step.end());
      pathSolutionIndex.erase(pathSolutionIndex.begin()+i+1,pathSolutionIndex.end());
      pathSolutionCoefficient.erase(pathSolutionCoefficient.begin()+i+1,pathSolutionCoefficient.end());
      break;
    }

    //reinitialize the solver for next lambda
    betaTemp = lasso.beta();
    betaTemp += 10*thresholdC;
    lasso.initializeBeta(betaTemp);

  }
  //t2=clock();


  return List::create(Named("variable")=wrap(pathSolutionIndex),Named("coefficient")=wrap(pathSolutionCoefficient),Named("lambda")=wrap(lambdaC)
                      , Named("mu")=wrap(mu), Named("step")=wrap(step));

}


RcppExport SEXP EMfusedLasso(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP lambda1, SEXP lambda2, SEXP intercept, SEXP maxStep, SEXP burn, SEXP eps, SEXP eps0, SEXP epsCG)
{
  //double t1,t2;
  //t1=clock();

  //convert parameters
  int p(as<int>(nbVar)), n(as<int>(nbIndiv)), maxStepC(as<int>(maxStep)), burnC(as<int>(burn));
  Real epsC(as<Real>(eps)), eps0C(as<Real>(eps0)), epsCGC(as<Real>(epsCG)), lambda1C(as<Real>(lambda1)), lambda2C(as<Real>(lambda2));
  bool interceptC=as<bool>(intercept);

  CArrayXX x(n,p);
  convertToArray(data,x);
  CVectorX y(n);
  convertToVector(response,nbIndiv,y);


  //center the data for intercept
  STK::CVectorX muX(p);
  STK::Real mu = 0;

  if(interceptC)
  {
    mu = y.sum()/n;
    y -= mu;

    for(int j = 1; j <= p; j++)
    {
      muX[j] = x.col(j).sum()/n;
      for(int i=1; i <= n; i++)
        x(i,j) -= muX[j];
    }
  }

  //create EM
  EM algo(maxStepC,burnC,epsC);

  //create fused lasso
  FusedLasso fusedlasso(&x,&y, lambda1C, lambda2C, eps0C, epsCGC);

  //run
  algo.run(&fusedlasso);

  std::vector<STK::Real> beta(p);
  for(int i = 1; i <= p; i++)
    beta[i-1]=fusedlasso.beta(i);

  //number of step made by the em
  int step = algo.step();

  return List::create(Named("coefficient")=wrap(beta), Named("mu")=wrap(mu), Named("step")=wrap(step));

}


RcppExport SEXP EMlogisticLasso(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP lambda, SEXP intercept, SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG)
{
  //double t1,t2;
  //t1=clock();

  //convert parameters
  int p(as<int>(nbVar)), n(as<int>(nbIndiv)), maxStepC(as<int>(maxStep)), burnC(as<int>(burn));
  Real epsC(as<Real>(eps)), thresholdC(as<Real>(threshold)), epsCGC(as<Real>(epsCG));
  bool interceptC=as<bool>(intercept);

  CArrayXX x(n,p);
  convertToArray(data,x);
  CVectorX y(n);
  convertToVector(response,nbIndiv,y);

  vector<STK::Real> lambdaC = as<vector<STK::Real> >(lambda);

  if(interceptC)
    interceptC=true;

  //center the data for intercept
//  STK::CVectorX muX(p);
//  STK::Real mu = 0;
//
//  if(interceptC)
//  {
//    mu = y.sum()/n;
//    y -= mu;
//
//    for(int j = 1; j <= p; j++)
//    {
//      muX[j] = x.col(j).sum()/n;
//      for(int i=1; i <= n; i++)
//        x(i,j) -= muX[j];
//    }
//  }


  //if lambdaC[0]=-1, we have to generate the lambda sequence with the same way as the glmnet package
  bool genLambda = false;
  if(lambdaC[0] == -1)
  {
    STK::CVectorX Xty(p);
    Xty=x.transpose() * y;

    genLambda = true;
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
  vector<vector<int> > pathSolutionIndex(lambdaC.size());
  vector<vector<STK::Real> > pathSolutionCoefficient(lambdaC.size());
  vector<int> step(lambdaC.size());
  vector<int> indexTemp;
  vector<STK::Real> coefficientTemp;

  //t2=clock();
  //create EM
  EM algo(maxStepC,burnC,epsC);
  LogisticLasso lasso( &x, &y, lambdaC[0], thresholdC, epsCGC);

  //run for all lambda
  STK::CVectorX betaTemp(p);
  for(int i = 0; i < (int) lambdaC.size(); i++)
  {
    //change the lambda
    lasso.setLambda(lambdaC[i]);

    //run algorithm
    Chrono::start();
    algo.run(&lasso);
    //t1=Chrono::elapsed();

    //stock the new values
    copySTKVectorInSTDVector(lasso.currentBeta(), coefficientTemp);
    pathSolutionCoefficient[i]=coefficientTemp;
    copySTKArray2DVectorInSTDVector(lasso.currentSet(), indexTemp);
    pathSolutionIndex[i]=indexTemp;
    step[i]=algo.step();

    //if there is 0 non-zeros coefficients, we stop at this values of lambda
    if( (coefficientTemp.size()==0))
    {
      lambdaC.erase(lambdaC.begin()+i+1,lambdaC.end());
      step.erase(step.begin()+i+1,step.end());
      pathSolutionIndex.erase(pathSolutionIndex.begin()+i+1,pathSolutionIndex.end());
      pathSolutionCoefficient.erase(pathSolutionCoefficient.begin()+i+1,pathSolutionCoefficient.end());
      break;
    }

    //reinitialize the solver for next lambda
    betaTemp = lasso.beta();
    betaTemp += 10*thresholdC;
    lasso.initializeBeta(betaTemp);

  }
  //t2=clock();


  return List::create(Named("variable")=wrap(pathSolutionIndex),Named("coefficient")=wrap(pathSolutionCoefficient),Named("lambda")=wrap(lambdaC)
                      , Named("mu")=wrap(0), Named("step")=wrap(step));

}


RcppExport SEXP EMlogisticFusedLasso(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP lambda1, SEXP lambda2, SEXP intercept, SEXP maxStep, SEXP burn, SEXP eps, SEXP eps0, SEXP epsCG)
{
  //double t1,t2;
  //t1=clock();

  //convert parameters
  int p(as<int>(nbVar)), n(as<int>(nbIndiv)), maxStepC(as<int>(maxStep)), burnC(as<int>(burn));
  Real epsC(as<Real>(eps)), eps0C(as<Real>(eps0)), epsCGC(as<Real>(epsCG)), lambda1C(as<Real>(lambda1)), lambda2C(as<Real>(lambda2));
  bool interceptC=as<bool>(intercept);

  CArrayXX x(n,p);
  convertToArray(data,x);
  CVectorX y(n);
  convertToVector(response,nbIndiv,y);
  
  if(interceptC)
    interceptC=true;

  //center the data for intercept
//  STK::CVectorX muX(p);
//  STK::Real mu = 0;
//
//  if(interceptC)
//  {
//    mu = y.sum()/n;
//    y -= mu;
//
//    for(int j = 1; j <= p; j++)
//    {
//      muX[j] = x.col(j).sum()/n;
//      for(int i=1; i <= n; i++)
//        x(i,j) -= muX[j];
//    }
//  }

  //create EM
  EM algo(maxStepC,burnC,epsC);

  //create fused lasso
  LogisticFusedLasso fusedlasso(&x,&y, lambda1C, lambda2C, eps0C, epsCGC);

  //run
  algo.run(&fusedlasso);

  std::vector<STK::Real> beta(p);
  for(int i = 1; i <= p; i++)
    beta[i-1]=fusedlasso.beta(i);

  //number of step made by the em
  int step = algo.step();

  return List::create(Named("coefficient")=wrap(beta), Named("mu")=wrap(0), Named("step")=wrap(step));

}


RcppExport SEXP EMCVLasso(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP lambda, SEXP nbFolds, SEXP intercept, SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG)
{
  //double t1,t2;
  //t1=clock();

  //convert parameters from R to C
  int p(as<int>(nbVar)), n(as<int>(nbIndiv)), maxStepC(as<int>(maxStep)), burnC(as<int>(burn)), nbFoldsC(as<int>(nbFolds));
  Real epsC(as<Real>(eps)), thresholdC(as<Real>(threshold)), epsCGC(as<Real>(epsCG));
  bool interceptC=as<bool>(intercept);

  CArrayXX x(n,p);
  convertToArray(data,x);
  CVectorX y(n);
  convertToVector(response,nbIndiv,y);
  vector<STK::Real> lambdaC = as<vector<STK::Real> >(lambda);

  //intercept in the model?
  if(interceptC)
  {
    STK::CVectorX muX(p);
    STK::Real mu = y.sum()/n;
    y -= mu;

    for(int j = 1; j <= p; j++)
    {
      muX[j] = x.col(j).sum()/n;
      for(int i=1; i <= n; i++)
        x(i,j) -= muX[j];
    }
  }

  //if lambdaC[0]=-1, we have to generate the lambda sequence with the same way as the glmnet package
  if(lambdaC[0] == -1)
  {
    //Xty
    STK::CVectorX Xty(p);
    Xty=x.transpose() * y;

    lambdaC.resize(100);
    int pos;
    //corrMax is the first lambda value in the lars sequence
    STK::Real corrMax =Xty.abs().maxElt(pos);
    STK::Real lambdaMinRatio, minLambda;
    if(n < p)
      lambdaMinRatio=0.01;
    else
      lambdaMinRatio=0.0001;

    minLambda = corrMax * lambdaMinRatio;
    STK::Real gapLambda = (log(corrMax)-log(minLambda))/99;
    lambdaC[99] = corrMax;
    for(int i = 99; i>0 ; i--)
      lambdaC[i-1] = exp(log(lambdaC[i]) - gapLambda);
  }

  Residuals measure;

  //create cv for lasso
  CVLasso<Lasso> lassocv;
  //set data
  lassocv.setX(x);
  lassocv.setY(y);
  //set cv parameters
  lassocv.setNbFolds(nbFoldsC);
  lassocv.setIndex(lambdaC);
  //set em parameters
  lassocv.setBurn(burnC);
  lassocv.setMaxStep(maxStepC);
  lassocv.setEps(epsC);
  //set CG parameter
  lassocv.setEpsCG(epsCGC);
  //set threshold for lassosolver
  lassocv.setThreshold(thresholdC);
  //set type of measure
  lassocv.setTypeMeasure(&measure);
  //initialize the class
  lassocv.initialize();

  //run cv
  lassocv.run2();

  //find the position of the lambda with the smallest cv error
  STK::Real minCV;
  int pos;
  minCV=lassocv.cv().minElt(pos);


  //t1=clock();
  //conversion from STK to std
  vector<double> cv(lambdaC.size()),cvError(lambdaC.size());
  STK::CVectorX stkCv(lambdaC.size()),stkCvError(lambdaC.size());
  stkCv=lassocv.cv();
  stkCvError=lassocv.cvError();

  for(int i = 1; i <= (int) lambdaC.size();i++)
  {
    cv[i-1]=stkCv[i];
    cvError[i-1]=stkCvError[i];
  }
  //t2=clock();


  return List::create(Named("lambda")=wrap(lambdaC), Named("cv")=wrap(cv), Named("cvError")=wrap(cvError), Named("minCV")=wrap(minCV),Named("lambda.optimal")=wrap(lambdaC[pos-1]));
}



RcppExport SEXP EMCVFusedLasso1D(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP lambda1, SEXP lambda2, SEXP optimL1, SEXP nbFolds, SEXP intercept, SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG)
{
  //double t1,t2;
  //t1=clock();

  //convert parameters
  int p(as<int>(nbVar)), n(as<int>(nbIndiv)), maxStepC(as<int>(maxStep)), burnC(as<int>(burn)), nbFoldsC(as<int>(nbFolds));
  Real epsC(as<Real>(eps)), thresholdC(as<Real>(threshold)), epsCGC(as<Real>(epsCG));
  bool interceptC=as<bool>(intercept), optimL1C=as<bool>(optimL1);

  CArrayXX x(n,p);
  convertToArray(data,x);
  CVectorX y(n);
  convertToVector(response,nbIndiv,y);

  vector<STK::Real> lambda1C = as<vector<STK::Real> >(lambda1);
  vector<STK::Real> lambda2C = as<vector<STK::Real> >(lambda2);

  //center the data
  if(interceptC)
  {
    STK::CVectorX muX(p);
    STK::Real mu = y.sum()/n;
    y -= mu;

    for(int j = 1; j <= p; j++)
    {
      muX[j] = x.col(j).sum()/n;
      for(int i=1; i <= n; i++)
        x(i,j) -= muX[j];
    }
  }

  //if lambda1 has to be optimized, we can generate the lambda1 sequence
  //if lambdaC[0]=-1, we have to generate the lambda sequence with the same way as the glmnet package
  if(optimL1C && (lambda1C[0] == -1) )
  {
    //Xty
    STK::CVectorX Xty(p);
    Xty=x.transpose() * y;

    lambda1C.resize(100);
    int pos;
    //corrMax is the first lambda value in the lars sequence
    STK::Real corrMax =Xty.abs().maxElt(pos);
    STK::Real lambdaMinRatio, minLambda;
    if(n < p)
      lambdaMinRatio=0.01;
    else
      lambdaMinRatio=0.0001;

    minLambda = corrMax * lambdaMinRatio;
    STK::Real gapLambda = (log(corrMax)-log(minLambda))/99;
    lambda1C[99] = corrMax;
    for(int i = 99; i>0 ; i--)
      lambda1C[i-1] = exp(log(lambda1C[i]) - gapLambda);
  }


  // lambda sequence to optimize
  vector<STK::Real> lambdaC;
  if(optimL1C)
    lambdaC = lambda1C;
  else
    lambdaC = lambda2C;

  //create the 1Dcv for fused lasso
  CVFusedLasso1D<FusedLasso> fusedlassocv;
  Residuals measure;
  //set the data
  fusedlassocv.setX(x);
  fusedlassocv.setY(y);
  //set the cv parameter
  fusedlassocv.setNbFolds(nbFoldsC);
  fusedlassocv.setIndex(lambdaC);
  if(optimL1C)
    fusedlassocv.setLambda(lambda2C[0]);//if we have to optimize lambda1, we set lambda2 as lambda parameter
  else
    fusedlassocv.setLambda(lambda1C[0]);
  //set EM parameters
  fusedlassocv.setBurn(burnC);
  fusedlassocv.setEps(epsC);
  fusedlassocv.setMaxStep(maxStepC);
  //set fusedlasso solver parameter
  fusedlassocv.setThreshold(thresholdC);
  //set CG parameter
  fusedlassocv.setEpsCG(epsCGC);
  fusedlassocv.setTypeMeasure(&measure);

  //initialize the created class
  fusedlassocv.initialize();

  //run the CV
  fusedlassocv.run2();

  //find the position of the lambda with the smallest cv error
  STK::Real minCV;
  int pos;
  minCV=fusedlassocv.cv().minElt(pos);

  //t1=clock();
  //conversion from STK to std vector
  int indexSize = lambdaC.size();

  vector<double> cv(indexSize),cvError(indexSize);
  STK::CVectorX stkCv, stkCvError;
  stkCv=fusedlassocv.cv();
  stkCvError=fusedlassocv.cvError();

  for(int i=1; i<=indexSize; i++)
  {
    cv[i-1]=stkCv[i];
    cvError[i-1]=stkCvError[i];
  }
  //t2=clock();

  return List::create(Named("lambda")=wrap(lambdaC), Named("cv")=wrap(cv), Named("cvError")=wrap(cvError), Named("minCV")=wrap(minCV),Named("lambda.optimal")=wrap(lambdaC[pos-1]));
}

RcppExport SEXP EMCVFusedLasso2D(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP lambda1, SEXP lambda2, SEXP nbFolds, SEXP intercept, SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG)
{
  //double t1,t2;
  //t1=clock();

  //convert parameters
  int p(as<int>(nbVar)), n(as<int>(nbIndiv)), maxStepC(as<int>(maxStep)), burnC(as<int>(burn)), nbFoldsC(as<int>(nbFolds));
  Real epsC(as<Real>(eps)), thresholdC(as<Real>(threshold)), epsCGC(as<Real>(epsCG));
  bool interceptC=as<bool>(intercept);

  CArrayXX x(n,p);
  convertToArray(data,x);
  CVectorX y(n);
  convertToVector(response,nbIndiv,y);

  vector<STK::Real> lambda1C = as<vector<STK::Real> >(lambda1);
  vector<STK::Real> lambda2C = as<vector<STK::Real> >(lambda2);

  //center the data
  if(interceptC)
  {
    STK::CVectorX muX(p);
    STK::Real mu = y.sum()/n;
    y -= mu;

    for(int j = 1; j <= p; j++)
    {
      muX[j] = x.col(j).sum()/n;
      for(int i=1; i <= n; i++)
        x(i,j) -= muX[j];
    }
  }

  //if lambdaC[0]=-1, we have to generate the lambda sequence with the same way as the glmnet package
  if(lambda1C[0] == -1)
  {
    //Xty
    STK::CVectorX Xty(p);
    Xty=x.transpose() * y;

    lambda1C.resize(100);
    int pos;
    //corrMax is the first lambda value in the lars sequence
    STK::Real corrMax =Xty.abs().maxElt(pos);
    STK::Real lambdaMinRatio, minLambda;
    if(n < p)
      lambdaMinRatio=0.01;
    else
      lambdaMinRatio=0.0001;

    minLambda = corrMax * lambdaMinRatio;
    STK::Real gapLambda = (log(corrMax)-log(minLambda))/99;
    lambda1C[99] = corrMax;
    for(int i = 99; i>0 ; i--)
      lambda1C[i-1] = exp(log(lambda1C[i]) - gapLambda);
  }

  //create the CV for fused lasso
  CVFusedLasso2D<FusedLasso> fusedlassocv;
  Residuals measure;
  //set the data
  fusedlassocv.setX(x);
  fusedlassocv.setY(y);
  //set CV parameters
  fusedlassocv.setNbFolds(nbFoldsC);
  fusedlassocv.setIndex(lambda1C);
  fusedlassocv.setIndexL2(lambda2C);
  //set EM parameters
  fusedlassocv.setEps(epsC);
  fusedlassocv.setMaxStep(maxStepC);
  fusedlassocv.setBurn(burnC);
  //set fusedsolver parameters
  fusedlassocv.setThreshold(thresholdC);
  //set CG parameter
  fusedlassocv.setEpsCG(epsCGC);
  fusedlassocv.setTypeMeasure(&measure);

  //initialize the class
  fusedlassocv.initialize();

  //run the algo
  fusedlassocv.run2();

  //find the position of the lambda with the smallest error
  STK::Real minCV;
  int pos;
  minCV=fusedlassocv.cv().minElt(pos);
  //convert the position in lambda1 and lambda2
  vector<STK::Real> lambdaMin(2);
  lambdaMin[0] = lambda1C[pos/lambda1C.size()];
  lambdaMin[1] = lambda2C[(pos%lambda2C.size())-1];

  //t1=clock();
  //conversion from stk to std
  int indexSize = lambda1C.size() * lambda2C.size();
  vector<double> cv(indexSize),cvError(indexSize);
  STK::CVectorX stkCv(indexSize),stkCvError(indexSize);
  stkCv=fusedlassocv.cv();
  stkCvError=fusedlassocv.cvError();

  for(int i=1; i<=indexSize; i++)
  {
    cv[i-1]=stkCv[i];
    cvError[i-1]=stkCvError[i];
  }
  //t2=clock();

  return List::create(Named("cv")=wrap(cv), Named("cvError")=wrap(cvError), Named("minCV")=wrap(minCV),Named("lambda.optimal")=wrap(lambdaMin));
}



// duplication pour cv logistic, trouver autre chose de mieux

RcppExport SEXP EMCVLogisticLasso(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP lambda, SEXP nbFolds, SEXP intercept, SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG)
{
  //double t1,t2;
  //t1=clock();

  //convert parameters from R to C
  int p(as<int>(nbVar)), n(as<int>(nbIndiv)), maxStepC(as<int>(maxStep)), burnC(as<int>(burn)), nbFoldsC(as<int>(nbFolds));
  Real epsC(as<Real>(eps)), thresholdC(as<Real>(threshold)), epsCGC(as<Real>(epsCG));
  bool interceptC=as<bool>(intercept);

  CArrayXX x(n,p);
  convertToArray(data,x);
  CVectorX y(n);
  convertToVector(response,nbIndiv,y);
  vector<STK::Real> lambdaC = as<vector<STK::Real> >(lambda);

  if(interceptC)
    interceptC=true;
    
  //intercept in the model?
//  if(interceptC)
//  {
//    STK::CVectorX muX(p);
//    STK::Real mu = y.sum()/n;
//    y -= mu;
//
//    for(int j = 1; j <= p; j++)
//    {
//      muX[j] = x.col(j).sum()/n;
//      for(int i=1; i <= n; i++)
//        x(i,j) -= muX[j];
//    }
//  }

  //if lambdaC[0]=-1, we have to generate the lambda sequence with the same way as the glmnet package
  if(lambdaC[0] == -1)
  {
    //Xty
    STK::CVectorX Xty(p);
    Xty=x.transpose() * y;

    lambdaC.resize(100);
    int pos;
    //corrMax is the first lambda value in the lars sequence
    STK::Real corrMax =Xty.abs().maxElt(pos);
    STK::Real lambdaMinRatio, minLambda;
    if(n < p)
      lambdaMinRatio=0.01;
    else
      lambdaMinRatio=0.0001;

    minLambda = corrMax * lambdaMinRatio;
    STK::Real gapLambda = (log(corrMax)-log(minLambda))/99;
    lambdaC[99] = corrMax;
    for(int i = 99; i>0 ; i--)
      lambdaC[i-1] = exp(log(lambdaC[i]) - gapLambda);
  }

  AUC measure;

  //create cv for lasso
  CVLasso<LogisticLasso> lassocv;
  //set data
  lassocv.setX(x);
  lassocv.setY(y);
  //set cv parameters
  lassocv.setNbFolds(nbFoldsC);
  lassocv.setIndex(lambdaC);
  //set em parameters
  lassocv.setBurn(burnC);
  lassocv.setMaxStep(maxStepC);
  lassocv.setEps(epsC);
  //set CG parameter
  lassocv.setEpsCG(epsCGC);
  //set threshold for lassosolver
  lassocv.setThreshold(thresholdC);
  //set type of measure
  lassocv.setTypeMeasure(&measure);
  //initialize the class
  lassocv.initialize();

  //run cv
  lassocv.run2();

  //find the position of the lambda with the smallest cv error
  STK::Real minCV;
  int pos;
  minCV=lassocv.cv().minElt(pos);


  //t1=clock();
  //conversion from STK to std
  vector<double> cv(lambdaC.size()),cvError(lambdaC.size());
  STK::CVectorX stkCv(lambdaC.size()),stkCvError(lambdaC.size());
  stkCv=lassocv.cv();
  stkCvError=lassocv.cvError();

  for(int i=1; i <= (int) lambdaC.size(); i++)
  {
    cv[i-1]=stkCv[i];
    cvError[i-1]=stkCvError[i];
  }
  //t2=clock();


  return List::create(Named("lambda")=wrap(lambdaC), Named("cv")=wrap(cv), Named("cvError")=wrap(cvError), Named("minCV")=wrap(minCV),Named("lambda.optimal")=wrap(lambdaC[pos-1]));
}



RcppExport SEXP EMCVLogisticFusedLasso1D(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP lambda1, SEXP lambda2, SEXP optimL1, SEXP nbFolds, SEXP intercept, SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG)
{
  //double t1,t2;
  //t1=clock();

  //convert parameters
  int p(as<int>(nbVar)), n(as<int>(nbIndiv)), maxStepC(as<int>(maxStep)), burnC(as<int>(burn)), nbFoldsC(as<int>(nbFolds));
  Real epsC(as<Real>(eps)), thresholdC(as<Real>(threshold)), epsCGC(as<Real>(epsCG));
  bool interceptC=as<bool>(intercept), optimL1C=as<bool>(optimL1);

  CArrayXX x(n,p);
  convertToArray(data,x);
  CVectorX y(n);
  convertToVector(response,nbIndiv,y);

  vector<STK::Real> lambda1C = as<vector<STK::Real> >(lambda1);
  vector<STK::Real> lambda2C = as<vector<STK::Real> >(lambda2);

  if(interceptC)
    interceptC=true;
    
  //center the data
//  if(interceptC)
//  {
//    STK::CVectorX muX(p);
//    STK::Real mu = y.sum()/n;
//    y -= mu;
//
//    for(int j = 1; j <= p; j++)
//    {
//      muX[j] = x.col(j).sum()/n;
//      for(int i=1; i <= n; i++)
//        x(i,j) -= muX[j];
//    }
//  }

  //if lambda1 has to be optimized, we can generate the lambda1 sequence
  //if lambdaC[0]=-1, we have to generate the lambda sequence with the same way as the glmnet package
  if(optimL1C && (lambda1C[0] == -1) )
  {
    //Xty
    STK::CVectorX Xty(p);
    Xty=x.transpose() * y;

    lambda1C.resize(100);
    int pos;
    //corrMax is the first lambda value in the lars sequence
    STK::Real corrMax =Xty.abs().maxElt(pos);
    STK::Real lambdaMinRatio, minLambda;
    if(n < p)
      lambdaMinRatio=0.01;
    else
      lambdaMinRatio=0.0001;

    minLambda = corrMax * lambdaMinRatio;
    STK::Real gapLambda = (log(corrMax)-log(minLambda))/99;
    lambda1C[99] = corrMax;
    for(int i = 99; i>0 ; i--)
      lambda1C[i-1] = exp(log(lambda1C[i]) - gapLambda);
  }


  // lambda sequence to optimize
  vector<STK::Real> lambdaC;
  if(optimL1C)
    lambdaC = lambda1C;
  else
    lambdaC = lambda2C;

  //create the 1Dcv for fused lasso
  CVFusedLasso1D<LogisticFusedLasso> fusedlassocv;
  AUC measure;
  //set the data
  fusedlassocv.setX(x);
  fusedlassocv.setY(y);
  //set the cv parameter
  fusedlassocv.setNbFolds(nbFoldsC);
  fusedlassocv.setIndex(lambdaC);
  if(optimL1C)
    fusedlassocv.setLambda(lambda2C[0]);//if we have to optimize lambda1, we set lambda2 as lambda parameter
  else
    fusedlassocv.setLambda(lambda1C[0]);
  //set EM parameters
  fusedlassocv.setBurn(burnC);
  fusedlassocv.setEps(epsC);
  fusedlassocv.setMaxStep(maxStepC);
  //set fusedlasso solver parameter
  fusedlassocv.setThreshold(thresholdC);
  fusedlassocv.setTypeMeasure(&measure);
  //set CG parameter
  fusedlassocv.setEpsCG(epsCGC);

  //initialize the created class
  fusedlassocv.initialize();

  //run the CV
  fusedlassocv.run2();

  //find the position of the lambda with the smallest cv error
  STK::Real minCV;
  int pos;
  minCV=fusedlassocv.cv().minElt(pos);

  //t1=clock();
  //conversion from STK to std vector
  int indexSize = lambdaC.size();

  vector<double> cv(indexSize),cvError(indexSize);
  STK::CVectorX stkCv, stkCvError;
  stkCv=fusedlassocv.cv();
  stkCvError=fusedlassocv.cvError();

  for(int i=1; i<=indexSize; i++)
  {
    cv[i-1]=stkCv[i];
    cvError[i-1]=stkCvError[i];
  }
  //t2=clock();

  return List::create(Named("lambda")=wrap(lambdaC), Named("cv")=wrap(cv), Named("cvError")=wrap(cvError), Named("minCV")=wrap(minCV),Named("lambda.optimal")=wrap(lambdaC[pos-1]));
}

RcppExport SEXP EMCVLogisticFusedLasso2D(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP lambda1, SEXP lambda2, SEXP nbFolds, SEXP intercept, SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG)
{
  //double t1,t2;
  //t1=clock();

  //convert parameters
  int p(as<int>(nbVar)), n(as<int>(nbIndiv)), maxStepC(as<int>(maxStep)), burnC(as<int>(burn)), nbFoldsC(as<int>(nbFolds));
  Real epsC(as<Real>(eps)), thresholdC(as<Real>(threshold)), epsCGC(as<Real>(epsCG));
  bool interceptC=as<bool>(intercept);

  CArrayXX x(n,p);
  convertToArray(data,x);
  CVectorX y(n);
  convertToVector(response,nbIndiv,y);

  vector<STK::Real> lambda1C = as<vector<STK::Real> >(lambda1);
  vector<STK::Real> lambda2C = as<vector<STK::Real> >(lambda2);

  if(interceptC)
    interceptC=true;
    
  //center the data
/*  if(interceptC)
  {
    STK::CVectorX muX(p);
    STK::Real mu = y.sum()/n;
    y -= mu;

    for(int j = 1; j <= p; j++)
    {
      muX[j] = x.col(j).sum()/n;
      for(int i=1; i <= n; i++)
        x(i,j) -= muX[j];
    }
  }*/

  //if lambdaC[0]=-1, we have to generate the lambda sequence with the same way as the glmnet package
  if(lambda1C[0] == -1)
  {
    //Xty
    STK::CVectorX Xty(p);
    Xty=x.transpose() * y;

    lambda1C.resize(100);
    int pos;
    //corrMax is the first lambda value in the lars sequence
    STK::Real corrMax =Xty.abs().maxElt(pos);
    STK::Real lambdaMinRatio, minLambda;
    if(n < p)
      lambdaMinRatio=0.01;
    else
      lambdaMinRatio=0.0001;

    minLambda = corrMax * lambdaMinRatio;
    STK::Real gapLambda = (log(corrMax)-log(minLambda))/99;
    lambda1C[99] = corrMax;
    for(int i = 99; i>0 ; i--)
      lambda1C[i-1] = exp(log(lambda1C[i]) - gapLambda);
  }

  //create the CV for fused lasso
  CVFusedLasso2D<LogisticFusedLasso> fusedlassocv;
  AUC measure;
  //set the data
  fusedlassocv.setX(x);
  fusedlassocv.setY(y);
  //set CV parameters
  fusedlassocv.setNbFolds(nbFoldsC);
  fusedlassocv.setIndex(lambda1C);
  fusedlassocv.setIndexL2(lambda2C);
  //set EM parameters
  fusedlassocv.setEps(epsC);
  fusedlassocv.setMaxStep(maxStepC);
  fusedlassocv.setBurn(burnC);
  //set fusedsolver parameters
  fusedlassocv.setThreshold(thresholdC);
  //set CG parameter
  fusedlassocv.setEpsCG(epsCGC);
  fusedlassocv.setTypeMeasure(&measure);
  //initialize the class
  fusedlassocv.initialize();

  //run the algo
  fusedlassocv.run2();

  //find the position of the lambda with the smallest error
  STK::Real minCV;
  int pos;
  minCV=fusedlassocv.cv().minElt(pos);
  //convert the position in lambda1 and lambda2
  vector<STK::Real> lambdaMin(2);
  lambdaMin[0] = lambda1C[pos/lambda1C.size()];
  lambdaMin[1] = lambda2C[(pos%lambda2C.size())-1];

  //t1=clock();
  //conversion from stk to std
  int indexSize = lambda1C.size() * lambda2C.size();
  vector<double> cv(indexSize),cvError(indexSize);
  STK::CVectorX stkCv(indexSize),stkCvError(indexSize);
  stkCv=fusedlassocv.cv();
  stkCvError=fusedlassocv.cvError();

  for(int i=1; i<=indexSize; i++)
  {
    cv[i-1]=stkCv[i];
    cvError[i-1]=stkCvError[i];
  }
  //t2=clock();

  return List::create(Named("cv")=wrap(cv), Named("cvError")=wrap(cvError), Named("minCV")=wrap(minCV),Named("lambda.optimal")=wrap(lambdaMin));
}
