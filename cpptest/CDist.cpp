#include "CDist.h"


CGaussianDist::CGaussianDist()
{
  setInitParam();
}
CGaussianDist::CGaussianDist(const CGaussianDist& dist) : precision(dist.precision)
{
}
CGaussianDist::~CGaussianDist()
{
}
// Gaussian prior.
void CGaussianDist::setParam(const double val, const int index)
{
  assert(index>=0);
  assert(index<getNumParams());
  switch(index)
    {
    case 0:
      precision = val;
      break;
    otherwise:
      cerr << "No such parameter" << endl;
      exit(1);
    }

}
double CGaussianDist::getParam(const int index) const
{
  assert(index>=0);
  assert(index<getNumParams());
  switch(index)
    {
    case 0:
      return precision;
    otherwise:
      cerr << "No such parameter" << endl;
      exit(1);
    }

}
void CGaussianDist::setInitParam()
{
  setType("gaussian");
  setName("Gaussian prior");
  setParamName("precision", 0);
  setNumParams(1);
  precision = 1.0;
  addTransform(new CNegLogLogitTransform, 0);
}
double CGaussianDist::logProb(double x) const
{
  return -0.5*precision*x*x -0.5*(ndlutil::LOGTWOPI - log(precision));
}
double CGaussianDist::getGradInput(double x) const
{
  return -precision*x;
}

CGammaDist::CGammaDist()
{
  setInitParam();
}
CGammaDist::CGammaDist(const CGammaDist& dist) : a(dist.a), b(dist.b)
{
}
CGammaDist::~CGammaDist()
{
}
// Gamma prior.
void CGammaDist::setParam(const double val, const int index)
{
  assert(index>=0);
  assert(index<getNumParams());
  switch(index)
    {
    case 0:
      a = val;
      break;
    case 1:
      b = val;
      break;
    otherwise:
      cerr << "No such parameter" << endl;
      exit(1);
    }

}
double CGammaDist::getParam(const int index) const
{
  assert(index>=0);
  assert(index<getNumParams());
  switch(index)
    {
    case 0:
      return a;
    case 1:
      return b;
    otherwise:
      cerr << "No such parameter" << endl;
      exit(1);
    }

}
void CGammaDist::setInitParam()
{
  setType("gamma");
  setName("gamma prior");
  setNumParams(2);  
  setParamName("a", 0);
  a=1e-6;
  setParamName("b", 1);
  b=1e-6;
  addTransform(new CNegLogLogitTransform, 0);
  addTransform(new CNegLogLogitTransform, 1);
}
double CGammaDist::logProb(double x) const
{
  
  return a*log(b) - ndlutil::gammaln(a) + ndlutil::xlogy(a-1.0,x)-b*x;  
}
double CGammaDist::getGradInput(double x) const
{
  return (a-1.0)/x - b;
}


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
  for(int i=0; i<getNumDists(); i++)
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
  int numDists = mxGetNumberOfElements(distArray);
  string distType;
  vector<int> distIndex;
  int counter = 0;
  for(int i=0; i<numDists; i++)
    {
      distType=mxArrayExtractStringField(distArray, "type", i);
      distIndex=mxArrayExtractVectorIntField(distArray, "index", i);
      for(int j=0; j<distIndex.size(); j++)
	{
	  counter++;
	  if(distType=="gamma")
	    addDist(new CGammaDist, distIndex[j]-1);
	  else if(distType=="gaussian")
	    addDist(new CGaussianDist, distIndex[j]-1);
	  else
	    cerr << "Dist type " << distType << " is currently unknown."<< endl;
	}
    }    
}
