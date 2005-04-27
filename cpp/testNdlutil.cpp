#include "ndlutil.h"
#include "CMatrix.h"
#include <iostream>

int testFunctions();
int testFunction(string funcName);
int testConsts();

int main()
{
  cout << "eps is: " << ndlutil::EPS << endl;
  int fail = 0;
  fail+= testFunctions();
  fail+= testConsts();
  cout << "Total failures: " << fail << endl;
}

int testConsts()
{
  int fail=0;
  double diff=abs(ndlutil::LOGTWOPI-log(2*M_PI));
  if(diff<ndlutil::EPS)
    cout << "log(2*pi) passes." << endl;
  else 
    {
      cout << "FAILURE: log(2*pi)." << ": absolute difference " << diff << endl;
      fail++;
    }
  diff=abs(ndlutil::HALFLOGTWOPI-.5*log(2*M_PI));
  if(diff<ndlutil::EPS)
    cout << ".5*log(2*pi) passes." << endl;
  else 
    {
      cout << "FAILURE: .5*log(2*pi)." << ": absolute difference " << diff << endl;
      fail++;
    }
  diff=abs(ndlutil::HALFSQRTTWO-.5*sqrt(2.0));
  if(diff<ndlutil::EPS)
    cout << ".5*sqrt(2) passes." << endl;
  else 
    {
      cout << "FAILURE: .5*sqrt(2)." << ": absolute difference " << diff << endl;
      fail++;
    }
  diff=abs(ndlutil::SQRTTWOPI-sqrt(2*M_PI));
  if(diff<ndlutil::EPS)
    cout << "sqrt(2*pi) passes." << endl;
  else 
    {
      cout << "FAILURE: sqrt(2*pi)." << ": absolute difference " << diff  << endl;
      fail++;
    }
  return fail;
}
//y=ndlutil::lnCumGaussSum(double u1, double u2, double w1, double w2);

int testFunctions()
{
  int fail = 0;
  fail+=testFunction("ngaussian");
  fail+=testFunction("cumGaussian");
  fail+=testFunction("invCumGaussian");
  fail+=testFunction("gradLnCumGaussian");
  fail+=testFunction("lnCumGaussian");
  fail+=testFunction("sigmoid");
  fail+=testFunction("invSigmoid");
  fail+=testFunction("erfcinv");
  fail+=testFunction("gamma");
  fail+=testFunction("gammaln");
  fail+=testFunction("digamma");
  return fail;
}
int testFunction(string funcName)
{
  int fail = 0;
  string fileName = funcName+"Test.mat";
  CMatrix xMat(1, 1);
  xMat.readMatlabFile(fileName, "x");
  CMatrix yMat(1, 1);
  yMat.readMatlabFile(fileName, "y");
  double x=xMat.getVal(0);
  double trueY=yMat.getVal(0);
  double y;
  if(funcName=="ngaussian")
    y=ndlutil::ngaussian(x);
  else if(funcName=="cumGaussian")    
    y=ndlutil::cumGaussian(x);
  else if(funcName=="gradLnCumGaussian")
    y=ndlutil::gradLnCumGaussian(x);
  else if(funcName=="invCumGaussian")
    y=ndlutil::invCumGaussian(x);
  else if(funcName=="lnCumGaussian")
    y=ndlutil::lnCumGaussian(x);
  else if(funcName=="sigmoid")
    y=ndlutil::sigmoid(x);
  else if(funcName=="invSigmoid")
    y=ndlutil::invSigmoid(x);
  else if(funcName=="erfcinv")
    y=ndlutil::erfcinv(x);
  else if(funcName=="gamma")
    y=ndlutil::gamma(x);
  else if(funcName=="gammaln")
    y=ndlutil::gammaln(x);
  else if(funcName=="digamma")
    y=ndlutil::digamma(x);
  double diff = abs(y-trueY);
  if(diff<=ndlutil::EPS)
    cout << funcName << " test passed." << endl;
  else
    {
      cout << "FAILURE: " << funcName << ": absolute difference " << diff << endl;
      fail++;
    }
  return fail;
}
  
