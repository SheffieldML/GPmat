#include "CTransform.h"


CNegLogLogitTransform::CNegLogLogitTransform()
{
  transform = 1;
  setType("negLogLogit");
}
  
double CNegLogLogitTransform::atox(double a) const
{
  double x;
  if(a<-limVal)
    x = exp(-limVal);
  else if(a<limVal)
    x=log(1+exp(a));
  assert(!isnan(x));
  return x;
}
double CNegLogLogitTransform::xtoa(double x) const
{
  double a=x;
  assert(a>-limVal);
  if(a<limVal)
    a=log(exp(x)-1);
  assert(!isnan(a));
  return a;
}
double CNegLogLogitTransform::gradfact(double x) const
{
  double g;
  assert(x>-limVal);
  if(x<limVal)
    g=(exp(x)-1)/exp(x);
  else
    g=1.0;
  return g;
}

CSigmoidTransform::CSigmoidTransform()
{
  transform = 1;
  setType("sigmoid");
}

double CSigmoidTransform::atox(double a) const
{
  double x;
  if(a<-limVal)
    return ndlutil::EPS;
  if(a<limVal)
    return ndlutil::sigmoid(a);
  else
    return 1-ndlutil::EPS;
}
double CSigmoidTransform::xtoa(double x) const
{
  return ndlutil::invSigmoid(x);
}
double CSigmoidTransform::gradfact(double x) const
{
  return x*(1-x);
}


mxArray* CParamTransforms::toMxArray() const
{
  int dims[1];
  // transforms field.
  const char *transFieldNames[] = {"index", "type"};
  dims[0]=getNumTransforms();
  mxArray* transformsArray = mxCreateStructArray(1, dims, 2, transFieldNames);
  CMatrix ind(1, 1);
  const char *compType[1];
  string trType;
  for(int i=0; i<getNumTransforms(); i++)
    {
      ind.setVals((double)(getTransformIndex(i)+1));
      mxSetField(transformsArray, i, "index", ind.toMxArray());
      trType = getTransformType(i);
      compType[0] = trType.c_str();
      mxSetField(transformsArray, i, "type", 
		 mxCreateCharMatrixFromStrings(1, compType));
    }
  return transformsArray;
}
void CParamTransforms::fromMxArray(const mxArray* transformArray)
{
  int numTransforms = mxGetNumberOfElements(transformArray);
  string transformType;
  vector<int> transformIndex;
  int counter = 0;
  for(int i=0; i<numTransforms; i++)
    {
      transformType=mxArrayExtractStringField(transformArray, "type", i);
      transformIndex=mxArrayExtractVectorIntField(transformArray, "index", i);
      for(int j=0; j<transformIndex.size(); j++)
	{
	  counter++;
	  if(transformType=="negLogLogit")
	    addTransform(new CNegLogLogitTransform, transformIndex[j]-1);
	  else if(transformType=="sigmoid")
	    addTransform(new CSigmoidTransform, transformIndex[j]-1);
	  else
	    cerr << "Transform type " << transformType << " is currently unknown."<< endl;
	}
    }    
}

