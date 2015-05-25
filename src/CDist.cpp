#include "CDist.h"


void CDist::writeParamsToStream(ostream& out) const
{
  writeToStream(out, "numParams", getNumParams());
  CMatrix par(1, getNumParams());
  getParams(par);
  par.toStream(out);
}
void CDist::readParamsFromStream(istream& in)
{
  string tbaseType = getBaseTypeStream(in);
  if(tbaseType != getBaseType())
    throw ndlexceptions::StreamFormatError("baseType", "Error mismatch between saved base type, " + tbaseType + ", and Class base type, " + getType() + ".");
  string ttype = getTypeStream(in);
  if(ttype != getType())
    throw ndlexceptions::StreamFormatError("type", "Error mismatch between saved type, " + ttype + ", and Class type, " + getType() + ".");
  unsigned int nPars = readIntFromStream(in, "numParams");
  CMatrix par(1, nPars);
  par.fromStream(in);
  if(nPars==getNumParams())
    setParams(par);
  else
    throw ndlexceptions::StreamFormatError("numParams", "Listed number of parameters does not match computed number of parameters.");
}
bool CDist::equals(const CDist& dist, double tol) const
{
  if(getType()!=dist.getType())
    return false;
  if(getNumParams()!=dist.getNumParams())
    return false;
  CMatrix params(1, getNumParams());
  getParams(params);
  CMatrix distParams(1, getNumParams());
  dist.getParams(distParams);
  if(!params.equals(distParams, tol))
    return false;
  return true;
}

#ifdef _NDLMATLAB
mxArray* CDist::toMxArray() const
{
  int dims[1];
  dims[0] = 1;
  const char *fieldNames[] = {"type", "transforms"};
  mxArray* matlabArray = mxCreateStructArray(1, dims, 2, fieldNames);
    
  // type field.
  const char *typeName[1];
  string ty=getType();
  typeName[0] = ty.c_str();
  mxSetField(matlabArray, 0, "type", 
	     mxCreateCharMatrixFromStrings(1, typeName));
  
  // transforms field.
  mxSetField(matlabArray, 0, "transforms", transformsToMxArray());
  
  // Class specific code.
  addParamToMxArray(matlabArray);
  return matlabArray;

}
void CDist::fromMxArray(const mxArray* matlabArray) 
{
  string mxType = mxArrayExtractStringField(matlabArray, "type");
  if(mxType!=getType())
  {
    throw ndlexceptions::FileReadError("Error mismatch between saved type, " + mxType + ", and Class type, " + getType() + ".");
  }
  mxArray* transformArray = mxArrayExtractMxArrayField(matlabArray, "transforms");
  // transforms field.
  transformsFromMxArray(transformArray);
  extractParamFromMxArray(matlabArray);
}
void CDist::extractParamFromMxArray(const mxArray* matlabArray)
{
  nParams = mxArrayExtractIntField(matlabArray, "nParams");
  string pName;
  for(unsigned int i=0; i<nParams; i++)
    {
      pName=getParamName(i);
      setParam(mxArrayExtractDoubleField(matlabArray, pName), i);  
    }
}
void CDist::addParamToMxArray(mxArray* matlabArray) const
{
  mxAddField(matlabArray, "nParams");
  mxSetField(matlabArray, 0, "nParams", convertMxArray((double)nParams));
  string pName;
  for(unsigned int i=0; i<nParams; i++)
    {      
      pName = getParamName(i);
      mxAddField(matlabArray, pName.c_str());      
      mxSetField(matlabArray, 0, pName.c_str(), convertMxArray(getParam(i))); 
    }
}
#endif /* _NDLMATLAB*/
CGaussianDist::CGaussianDist()
{
  _init();
  setInitParam();
}
CGaussianDist::CGaussianDist(const CGaussianDist& dist) : precision(dist.precision)
{
  _init();
}
CGaussianDist::~CGaussianDist()
{
}
// Gaussian prior.
void CGaussianDist::setParam(double val, unsigned int index)
{
  
  BOUNDCHECK(index<getNumParams());
  switch(index)
  {
  case 0:
    precision = val;
    break;
  default:
    throw ndlexceptions::Error("Requested parameter doesn't exist.");
  } 
}
double CGaussianDist::getParam(unsigned int index) const
{
  
  BOUNDCHECK(index<getNumParams());
  switch(index)
  {
  case 0:
    return precision;
  default:
    throw ndlexceptions::Error("Requested parameter doesn't exist.");
  }
}
void CGaussianDist::_init()
{
  setNumParams(1);
  setType("gaussian");
  setName("Gaussian prior");
  setParamName("precision", 0);
  addTransform(new CNegLogLogitTransform, 0);
}
void CGaussianDist::setInitParam()
{
  precision = 1.0;
}
double CGaussianDist::logProb(double x) const
{
  return -0.5*precision*x*x -0.5*(ndlutil::LOGTWOPI - log(precision));
}
double CGaussianDist::getGradInput(double x) const
{
  return -precision*x;
}

// Gamma prior.
CGammaDist::CGammaDist()
{
  _init();
  setInitParam();
}
CGammaDist::CGammaDist(const CGammaDist& dist) : a(dist.a), b(dist.b)
{
  _init();
}
CGammaDist::~CGammaDist()
{
}
void CGammaDist::setParam( double val,  unsigned int index)
{
  
  BOUNDCHECK(index<getNumParams());
  switch(index)
  {
  case 0:
    a = val;
    break;
  case 1:
    b = val;
    break;
  default:
    throw ndlexceptions::Error("Requested parameter doesn't exist.");
  }
}
double CGammaDist::getParam(unsigned int index) const
{
  
  BOUNDCHECK(index<getNumParams());
  switch(index)
  {
  case 0:
    return a;
  case 1:
    return b;
  default:
    throw ndlexceptions::Error("Requested parameter doesn't exist.");
  }
}
void CGammaDist::_init()
{
  setNumParams(2);  
  setType("gamma");
  setName("gamma prior");
  setParamName("a", 0);
  addTransform(new CNegLogLogitTransform, 0);
  setParamName("b", 1);
  addTransform(new CNegLogLogitTransform, 1);
}
void CGammaDist::setInitParam()
{
  a=1e-6;
  b=1e-6;
}
double CGammaDist::logProb(double x) const
{
  return a*log(b) - ndlutil::gammaln(a) + ndlutil::xlogy(a-1.0,x)-b*x;  
}
double CGammaDist::getGradInput(double x) const
{
  return (a-1.0)/x - b;
}

// Wang's unusual prior from the GPDM thesis.
CWangDist::CWangDist()
{
  _init();
  setInitParam();
}
CWangDist::CWangDist(const CWangDist& dist) : M(dist.M)
{
  _init();
}
CWangDist::~CWangDist()
{
}
void CWangDist::setParam( double val,  unsigned int index)
{
  
  BOUNDCHECK(index<getNumParams());
  switch(index)
  {
  case 0:
    M = val;
    break;
  default:
    throw ndlexceptions::Error("Requested parameter doesn't exist.");
  }
}
double CWangDist::getParam(unsigned int index) const
{
  
  BOUNDCHECK(index<getNumParams());
  switch(index)
  {
  case 0:
    return M;
  default:
    throw ndlexceptions::Error("Requested parameter doesn't exist.");
  }
}
void CWangDist::_init()
{
  setType("wang");
  setName("Wang's GPDM prior");
  setNumParams(1);  
  setParamName("M", 0);
  addTransform(new CNegLogLogitTransform, 0);
}
void CWangDist::setInitParam()
{
  M=1;
}
double CWangDist::logProb(double x) const
{
  return -M*log(x);
}
double CWangDist::getGradInput(double x) const
{
  return -M/x;
}

#ifdef _NDLMATLAB
mxArray* CParamPriors::toMxArray() const
{
  int dims[1];
  // dists field.
  const char *transFieldNames[] = {"index", "type"};
  dims[0]=getNumDists();
  mxArray* distsArray = mxCreateStructArray(1, dims, 2, transFieldNames);
  CMatrix ind(1, 1);
  const char *compType[1];
  string trType;
  for(unsigned int i=0; i<getNumDists(); i++)
  {
    ind.setVals((double)(getDistIndex(i)+1));
    mxSetField(distsArray, i, "index", ind.toMxArray());
    trType = getDistType(i);
    compType[0] = trType.c_str();
    mxSetField(distsArray, i, "type", 
	       mxCreateCharMatrixFromStrings(1, compType));
  }
  return distsArray;
}
void CParamPriors::fromMxArray(const mxArray* distArray)
{
  unsigned int numDists = mxGetNumberOfElements(distArray);
  string distType;
  vector<unsigned int> distIndex;
  unsigned int counter = 0;
  for(unsigned int i=0; i<numDists; i++)
  {
    distType=mxArrayExtractStringField(distArray, "type", i);
    distIndex=mxArrayExtractVectorUintField(distArray, "index", i);
    for(unsigned int j=0; j<distIndex.size(); j++)
    {
      counter++;
      if(distType=="gamma")
	addDist(new CGammaDist, distIndex[j]-1);
      else if(distType=="wang")
	addDist(new CWangDist, distIndex[j]-1);
      else if(distType=="gaussian")
	addDist(new CGaussianDist, distIndex[j]-1);
      else
	throw ndlexceptions::NotImplementedError("Dist type " + distType + " is currently not implemented.");
    }
  }    
}

#endif

void writeDistToStream(const CDist& dist, ostream& out)
{
  dist.toStream(out);
}
CDist* readDistFromStream(istream& in)
{
  double ver = CStreamInterface::readVersionFromStream(in); 
  string tbaseType = CStreamInterface::getBaseTypeStream(in);
  if(tbaseType != "dist")
    throw ndlexceptions::StreamFormatError("baseType", "Error mismatch between saved base type, " + tbaseType + ", and Class base type, dist.");
  CDist* pdist;

  string type = CStreamInterface::getTypeStream(in);
  if(type=="gaussian")
    pdist = new CGaussianDist();
  else if(type=="gamma")
    pdist = new CGammaDist();
  else if(type=="wang")
    pdist = new CWangDist();
  else
    throw ndlexceptions::StreamFormatError("type", "Unknown distribution type " + type + ".");
  pdist->readParamsFromStream(in);  
  return pdist;
}
