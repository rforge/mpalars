#include "stkpp/projects/STKernel/include/STK_Real.h"
#include "larsR.h"
#include "lars/Lars.h"
#include "fusion/Fusion.h"

#include <iostream>

using namespace Rcpp;
using namespace STK;
using namespace std;
using namespace MPA;

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

RcppExport SEXP lars(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP maxStep, SEXP eps, SEXP verbose)
{
  double t1,t2;
  t1=clock();
  //convert parameters
  int p(as<int>(nbVar)), n(as<int>(nbIndiv)), maxStepC(as<int>(maxStep));
  bool print=as<bool>(verbose);
  Real epsC(as<Real>(eps));

  CArrayXX x(n,p);
  convertToArray(data,x);
  CVectorX y(n);
  convertToVector(response,nbIndiv,y);

  t2=clock();
  if(print)
    cout<<"Temps conversion des données:"<<(double) (t2-t1)/CLOCKS_PER_SEC<<"s"<<endl;

  //run algorithm
  t1=clock();
  Lars lars(x,y,maxStepC,epsC,print);
  lars.run();
  t2=clock();

  if(print)
    cout<<"Temps lars:"<<(double) (t2-t1)/CLOCKS_PER_SEC<<"s"<<endl;

  //extract and convert results
  t1=clock();
  int step=lars.step();
  vector<double> lambda(step+1);
  vector<vector<int> > varIdx(step+1);
  vector<vector<double> > varCoeff(step+1);
  vector<int> evoIdxDrop(step,0);
  vector<int> evoIdxAdd(step,0);

  lambda[0]=0;
  for(int i=1;i<=step;i++)
  {
    varIdx[i].resize(lars.path(i).size());
    varCoeff[i].resize(lars.path(i).size());
    for(int j=1;j<=lars.path(i).size();j++)
    {
      varCoeff[i][j-1]=lars.coefficient(i,j);
      varIdx[i][j-1]=lars.varIdx(i,j);
    }
    lambda[i]=lars.lambda(i);
    if(lars.evolution()[i-1].first!=0)
        evoIdxAdd[i-1]=lars.evolution()[i-1].first;
    if(lars.evolution()[i-1].second!=0)
        evoIdxDrop[i-1]=lars.evolution()[i-1].second;

  }
  vector<int> ignored;
  ignored.reserve(p);
  for(int i=1;i<=lars.toIgnore().size();i++)
	if(lars.toIgnore()[i])
		ignored.push_back(i);


  t2=clock();
  if(print)
    cout<<"Temps extraction des données:"<<(double) (t2-t1)/CLOCKS_PER_SEC<<"s"<<endl;

  return List::create(Named("lambda")=wrap(lambda),Named("varIdx")=wrap(varIdx),Named("varCoeff")=wrap(varCoeff),
                      Named("evoDropIdx")=wrap(evoIdxDrop), Named("evoAddIdx")=wrap(evoIdxAdd),Named("step")=wrap(step),Named("mu")=wrap(lars.mu()),Named("ignored")=wrap(ignored));

}

RcppExport SEXP fusion(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP maxStep, SEXP eps, SEXP verbose)
{
  double t1,t2;
  t1=clock();
  //convert parameters
  int p(as<int>(nbVar)), n(as<int>(nbIndiv)), maxStepC(as<int>(maxStep));
  bool print=as<bool>(verbose);
  Real epsC(as<Real>(eps));

  CArrayXX x(n,p);
  convertToArray(data,x);
  CVectorX y(n,p);
  convertToVector(response,nbIndiv,y);

  t2=clock();
  if(print)
    cout<<"Temps conversion des données:"<<(double) (t2-t1)/CLOCKS_PER_SEC<<"s"<<endl;

  //run algorithm
  t1=clock();
  Fusion fusion(x,y,maxStepC,epsC,print);
  fusion.run();
  t2=clock();

  if(print)
    cout<<"Temps fusion:"<<(double) (t2-t1)/CLOCKS_PER_SEC<<"s"<<endl;

  //extract and convert results
  t1=clock();
  int step=fusion.step();


  vector<double> lambda(step+1);
  vector<vector<int> > varIdx(step+1);
  vector<vector<double> > varCoeff(step+1);

  vector<int> evoIdxDrop(step,0);
  vector<int> evoIdxAdd(step,0);

  lambda[0]=0;
  for(int i=1;i<=step;i++)
  {
    varIdx[i].resize(fusion.path(i).size());
    varCoeff[i].resize(fusion.path(i).size());
    for(int j=1;j<=fusion.path(i).size();j++)
    {
      varCoeff[i][j-1]=fusion.coefficient(i,j);
      varIdx[i][j-1]=fusion.varIdx(i,j);
    }
    lambda[i]=fusion.lambda(i);
    if(fusion.evolution()[i-1].first!=0)
        evoIdxAdd[i-1]=fusion.evolution()[i-1].first;
    if(fusion.evolution()[i-1].second!=0)
        evoIdxDrop[i-1]=fusion.evolution()[i-1].second;

  }
  vector<int> ignored;
  ignored.reserve(p);
  for(int i=1;i<=fusion.toIgnore().size();i++)
	if(fusion.toIgnore()[i])
		ignored.push_back(i);
  t2=clock();
  if(print)
    cout<<"Temps extraction des données:"<<(double) (t2-t1)/CLOCKS_PER_SEC<<"s"<<endl;

  return List::create(Named("lambda")=wrap(lambda), Named("varIdx")=wrap(varIdx),Named("varCoeff")=wrap(varCoeff),
                      Named("step")=wrap(step),Named("evoDropIdx")=wrap(evoIdxDrop), Named("evoAddIdx")=wrap(evoIdxAdd),Named("mu")=wrap(fusion.mu()),Named("ignored")=wrap(ignored));

}
