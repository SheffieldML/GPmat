#include "CDist.h"

CDist::CDist()
{
}
CDist::~CDist()
{
}
/*CDist::CDist(CDist& dist)
{
  transforms = dist.transforms;
}*/


// Gaussian prior.
void CGaussian::setParams(CMatrix& params)
{
  precision = params.getVals(0);
}
CMatrix CGaussian::getParams()
{
  CMatrix params(precision);
}
void CGaussian::setInitParam()
{
  nParams = 1;
  precision = 1;
  addTransform(new CNegLogLogitTransform, 1);
}
double CGaussian::logProb(CMatrix& X)
{
  double x2 = 0;
  int numElements = X.getNumElements();
  for(int i = 0; i < numElements; i++)
    x2 += X.getVals(i)*X.getVals(i);
  return -0.5*precision*x2 -0.5*numElements*(log(2*M_PI) - log(precision));
}
CMatrix CGaussian::gradParam(CMatrix& X)
{
  CMatrix G(X);
  G.multiply(-precision);
  return G;
}
