#include "stkpp/projects/STKernel/include/STK_Real.h"
#include "larsR.h"
#include "HDPenReg/lars/Lars.h"
#include "HDPenReg/lars/Cvlars.h"
#include "HDPenReg/lars/Fusion.h"

#include <iostream>

using namespace Rcpp;
using namespace STK;
using namespace std;
using namespace HD;

void convertToVector(SEXP const& rVector, SEXP size, CVectorX &output)
{
  NumericVector data(rVector);
  int n=as<int>(size);

  for(int i=0;i<n;i++)
    output[i+1]=data[i];
}

void convertToArray(SEXP const& rMatrix, CArrayXX &output)
{
  NumericMatrix data(rMatrix);
  int n = data.nrow(), p=data.ncol();

  for(int i=1;i<=n;i++)
    for(int j=1;j<=p;j++)
      output(i,j)=data[(j-1)*n+i-1];
}

RcppExport SEXP lars(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP maxStep, SEXP intercept, SEXP eps)
{
  //double t1,t2;
  //t1=clock();
  //convert parameters
  int p(as<int>(nbVar)), n(as<int>(nbIndiv)), maxStepC(as<int>(maxStep));
  bool interceptC = as<bool>(intercept);
  Real epsC(as<Real>(eps));

  CArrayXX x(n,p);
  convertToArray(data,x);
  CVectorX y(n);
  convertToVector(response,nbIndiv,y);

  //t2=clock();

  //cout<<"Temps conversion des données:"<<(double) (t2-t1)/CLOCKS_PER_SEC<<"s"<<endl;

  //run algorithm
  //t1=clock();
  Lars lars(x,y,maxStepC,interceptC,epsC);
  lars.run();
  //t2=clock();

  //cout<<"Temps lars:"<<(double) (t2-t1)/CLOCKS_PER_SEC<<"s"<<endl;

  //extract and convert results
  //t1=clock();
  int step=lars.step();
  vector<double> l1norm(step+1);
  vector<STK::Real> lambda=lars.lambda();
  vector<vector<int> > varIdx(step+1);
  vector<vector<double> > varCoeff(step+1);
  vector<vector<int> > evoIdxDrop(step);
  vector<vector<int> > evoIdxAdd(step);
  vector<double> muX(p);

  for(int i = 1; i<= p; i++)
    muX[i-1] = lars.muX(i);
    
  l1norm[0]=0;
  for(int i=1;i<=step;i++)
  {
    varIdx[i].resize(lars.path(i).size());
    varCoeff[i].resize(lars.path(i).size());
    for(int j=1;j<=lars.path(i).size();j++)
    {
      varCoeff[i][j-1]=lars.coefficient(i,j);
      varIdx[i][j-1]=lars.varIdx(i,j);
    }
    l1norm[i]=lars.l1norm(i);
    if(lars.evolution()[i-1].first.size()!=0)
        evoIdxAdd[i-1]=lars.evolution()[i-1].first;
    if(lars.evolution()[i-1].second.size()!=0)
        evoIdxDrop[i-1]=lars.evolution()[i-1].second;

  }
  vector<int> ignored;
  ignored.reserve(p);
  for(int i = 0; i < (int) lars.toIgnore().size(); i++)
  {
    if(lars.toIgnore()[i])
  	  ignored.push_back(i);
  }



  //t2=clock();
  //cout<<"Temps extraction des données:"<<(double) (t2-t1)/CLOCKS_PER_SEC<<"s"<<endl;

  return List::create(Named("l1norm")=wrap(l1norm), Named("lambda")=wrap(lambda), Named("varIdx")=wrap(varIdx), Named("varCoeff")=wrap(varCoeff),
                      Named("evoDropIdx")=wrap(evoIdxDrop), Named("evoAddIdx")=wrap(evoIdxAdd),Named("step")=wrap(step),Named("mu")=wrap(lars.mu()),
                      Named("ignored")=wrap(ignored),Named("error")=wrap(lars.msg_error()),Named("muX")=muX);

}

RcppExport SEXP fusion(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP maxStep, SEXP intercept, SEXP eps)
{
  //double t1,t2;
  //t1=clock();
  //convert parameters
  int p(as<int>(nbVar)), n(as<int>(nbIndiv)), maxStepC(as<int>(maxStep));
  bool interceptC = as<bool>(intercept);
  Real epsC(as<Real>(eps));

  CArrayXX x(n,p);
  convertToArray(data,x);
  CVectorX y(n,p);
  convertToVector(response,nbIndiv,y);

  //t2=clock();
  //cout<<"Temps conversion des données:"<<(double) (t2-t1)/CLOCKS_PER_SEC<<"s"<<endl;

  //run algorithm
  //t1=clock();
  Fusion fusion(x,y,maxStepC,interceptC,epsC);
  fusion.run();
  //t2=clock();

  //cout<<"Temps fusion:"<<(double) (t2-t1)/CLOCKS_PER_SEC<<"s"<<endl;

  //extract and convert results
  //t1=clock();
  int step=fusion.step();


  vector<double> l1norm(step+1);
  vector<vector<int> > varIdx(step+1);
  vector<vector<double> > varCoeff(step+1);
  vector<STK::Real> lambda=fusion.lambda();
  vector<vector<int> > evoIdxDrop(step);
  vector<vector<int> > evoIdxAdd(step);
  vector<double> muX(p);

  for(int i = 1; i<= p; i++)
    muX[i-1] = fusion.muX(i);
  l1norm[0]=0;
  for(int i = 1; i <= step; i++)
  {
    varIdx[i].resize(fusion.path(i).size());
    varCoeff[i].resize(fusion.path(i).size());
    for(int j = 1; j <= fusion.path(i).size(); j++)
    {
      varCoeff[i][j-1]=fusion.coefficient(i,j);
      varIdx[i][j-1]=fusion.varIdx(i,j);
    }
    l1norm[i]=fusion.l1norm(i);
    if(fusion.evolution()[i-1].first.size()!=0)
        evoIdxAdd[i-1]=fusion.evolution()[i-1].first;
    if(fusion.evolution()[i-1].second.size()!=0)
        evoIdxDrop[i-1]=fusion.evolution()[i-1].second;

  }
  vector<int> ignored;
  ignored.reserve(p);
  for(int i = 0; i < (int) fusion.toIgnore().size(); i++)
  {
    if(fusion.toIgnore()[i])
		  ignored.push_back(i);
  }

  //t2=clock();
  //cout<<"Temps extraction des données:"<<(double) (t2-t1)/CLOCKS_PER_SEC<<"s"<<endl;

  return List::create(Named("l1norm")=wrap(l1norm), Named("lambda")=wrap(lambda), Named("varIdx")=wrap(varIdx),Named("varCoeff")=wrap(varCoeff),
                      Named("step")=wrap(step),Named("evoDropIdx")=wrap(evoIdxDrop), Named("evoAddIdx")=wrap(evoIdxAdd),Named("mu")=wrap(fusion.mu()),
                      Named("ignored")=wrap(ignored),Named("error")=wrap(fusion.msg_error()),Named("muX")=muX);

}


RcppExport SEXP cvlars(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP maxStep, SEXP intercept, SEXP eps, SEXP nbFold, SEXP partition, SEXP index)
{
  //double t1,t2;
  //t1=clock();
  //convert parameters
  int p(as<int>(nbVar)), n(as<int>(nbIndiv)), maxStepC(as<int>(maxStep)), nbFoldC(as<int>(nbFold));
  vector<double> indexC=as<vector<double> >(index);
  vector<int> partitionC=as<vector<int> >(partition);
  bool interceptC = as<bool>(intercept);
  Real epsC(as<Real>(eps));

  CArrayXX x(n,p);
  convertToArray(data,x);
  CVectorX y(n);
  convertToVector(response,nbIndiv,y);

  //t2=clock();

  //run algorithm
  //t1=clock();
  Cvlars cvlars(x,y,nbFoldC,indexC,maxStepC,interceptC,epsC);
  if(partitionC[0]!=-1)
    cvlars.setPartition(partitionC);
    
#ifdef SUPPORT_OPENMP   
  cvlars.run2();
#else
  cvlars.run();
#endif
  //t2=clock();

  //extract and convert results
  //t1=clock();
  vector<double> cv(indexC.size()),cvError(indexC.size());
  STK::CVectorX stkCv(indexC.size()),stkCvError(indexC.size());
  stkCv=cvlars.cv();
  stkCvError=cvlars.cvError();

  for(int i = 1; i <= (int) indexC.size(); i++)
  {
    cv[i-1]=stkCv[i];
    cvError[i-1]=stkCvError[i];
  }
  //t2=clock();


  return List::create(Named("cv")=wrap(cv),Named("cvError")=wrap(cvError));

}
