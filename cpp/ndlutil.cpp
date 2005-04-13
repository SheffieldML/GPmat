#include "ndlutil.h"

double ngaussian(double x)
{
  x *= x;
  x = exp(-.5*x);
  x = x/sqrt(2.0*M_PI);
  return x;
}
double cumGaussian(double x)
{
  x *= sqrt(2.0)/2;
  x = erf(x);
  x++;
  x*=0.5;
  return x;
}
double gradLnCumGaussian(double x)
{
  if(x>0)
    x = ngaussian(x)/cumGaussian(x);
  else
    x = 1/(sqrt(2.0*M_PI)*0.5*derfcx_(-sqrt(2.0)/2.0*x));
  return x;
}
double lnCumGaussian(double x)
{
  if(x<0)
    x = -.5*x*x + log(0.5) + log(derfcx_(-sqrt(2.0)/2.0*x));
  else
    x = log(cumGaussian(x));
  return x;
}
