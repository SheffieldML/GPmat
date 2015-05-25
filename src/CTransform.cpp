#include "CTransform.h"

const double CTransform::eps = 1e-16;
CTransform* CTransform::defaultPositive()
{
  return new CExpTransform();
}
CTransform* CTransform::defaultZeroOne()
{
  return new CSigmoidTransform();
}
CTransform* CTransform::getNewTransformPointer(const string transformType)
{
  if(transformType=="negLogLogit")
    return new CNegLogLogitTransform();
  else if(transformType=="sigmoid")
    return new CSigmoidTransform();
  else if(transformType=="exp")
    return new CExpTransform();
  else if(transformType=="linear")
    return new CLinearTransform();
  else
    throw ndlexceptions::Error("Transform type " + transformType + " is currently unknown.");
}
CExpTransform::CExpTransform()
{
  transform = true;
  setType("exp");
}
  
double CExpTransform::atox(double a) const
{
  double x=0;
  if(a<-limVal)
    x = exp(-limVal);
  else if(a<limVal)
    x=exp(a);
  else
    x = exp(limVal);
  SANITYCHECK(!isnan(x));
  return x;
}
double CExpTransform::xtoa(double x) const
{
  SANITYCHECK(x>0);
  double a=log(x);
  SANITYCHECK(!isnan(a));
  return a;
}
double CExpTransform::gradfact(double x) const
{
  return x;
}

CNegLogLogitTransform::CNegLogLogitTransform()
{
  transform = 1;
  setType("negLogLogit");
}
  
double CNegLogLogitTransform::atox(double a) const
{
  double x=0;
  if(a<-limVal)
    x = exp(-limVal);
  else if(a<limVal)
    x=log(1+exp(a));
  SANITYCHECK(!isnan(x));
  return x;
}
double CNegLogLogitTransform::xtoa(double x) const
{
  double a=x;
  SANITYCHECK(a>-limVal);
  if(a<limVal)
    a=log(exp(x)-1);
  SANITYCHECK(!isnan(a));
  return a;
}
double CNegLogLogitTransform::gradfact(double x) const
{
  double g;
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

void CParamTransforms::writeParamsToStream(ostream& out) const
{
  writeToStream(out, "numTransforms", getNumTransforms());
  for(unsigned int i=0; i<getNumTransforms(); i++)
  {
    writeToStream(out, "type", getTransformType(i));
    writeToStream(out, "index", getTransformIndex(i));
  }
}
void CParamTransforms::readParamsFromStream(istream& in)
{
  transforms.clear();
  transIndex.clear();
  int numTrans = readIntFromStream(in, "numTransforms");
  for(int i=0; i<numTrans; i++)
  {
    string trtype = readStringFromStream(in, "type");
    int trIndex = readIntFromStream(in, "index");
    addTransform(CTransform::getNewTransformPointer(trtype), trIndex);
  }
}
void CParamTransforms::display(ostream& out) const
{
  out << "Parameter Transforms:" << endl;
  for(unsigned int i=0; i<getNumTransforms(); i++)
  {
    out << "Transform type: " << getTransformType(i) << endl;
    out << "Transform index: " << getTransformIndex(i) << endl;
  }
}
bool CParamTransforms::equals(CParamTransforms transforms) const
{
  if(getNumTransforms()!=transforms.getNumTransforms())
    return false;
  for(unsigned int i=0; i<getNumTransforms(); i++)
  {
    if(getTransformType(i)!=transforms.getTransformType(i))
      return false;
    if(getTransformIndex(i)!=transforms.getTransformIndex(i))
      return false;
  }
  return true;
}

#ifdef _NDLMATLAB
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
    mxSetField(transformsArray, i, "type", mxCreateCharMatrixFromStrings(1, compType));
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
      CTransform* trans = CTransform::getNewTransformPointer(transformType);
      addTransform(trans, transformIndex[j]-1);
    }
  }
}

#endif
