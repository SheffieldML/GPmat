#include "CKern.h"

using namespace std;

void CKern::initialiseKern(const int inDim) 
{
  inputDim=inDim;
}
void CKern::initialiseKern(const CMatrix& X)
{
  inputDim=X.getCols();
}
// Copy constructor
void CKern::initialiseKern(const CKern& kern)
{
  // TODO I think this is messed up ...copy vectors and priors
  // transforms=kern.transforms;
  // TODO  priors=kern.priors;
  inputDim=kern.inputDim;
}
ostream& CKern::display(ostream& os) const
{
  os << getName() << " kernel:" << endl;
  for(int i=0; i<nParams; i++)
    {
      os << getParamName(i) << ": " << getParam(i) << endl;
    }
}

mxArray* CKern::toMxArray() const
{
  int dims[1];
  dims[0] = 1;
  const char *fieldNames[] = {"type", "transforms", "inputDimension"};
  mxArray* matlabArray = mxCreateStructArray(1, dims, 3, fieldNames);
    
  // type field.
  const char *typeName[1];
  string ty=getType();
  typeName[0] = ty.c_str();
  mxSetField(matlabArray, 0, "type", 
	     mxCreateCharMatrixFromStrings(1, typeName));
  
  // transforms field.
  mxSetField(matlabArray, 0, "transforms", transformsToMxArray());

  // inputDimension field.
  mxSetField(matlabArray, 0, "inputDimension", convertMxArray((double)inputDim));

  
  // whiteVariance field.
  double whiteVar = getWhite();
  if(whiteVar !=0.0)
    {
      mxAddField(matlabArray, "whiteVariance");      
      mxSetField(matlabArray, 0, "whiteVariance", convertMxArray(whiteVar)); 
    }
  
  // priors field.
  /// if priros exist need to add the field.

  // Kernel specific code.
  addParamToMxArray(matlabArray);
  return matlabArray;

}
void CKern::fromMxArray(const mxArray* matlabArray) 
{
  string mxType = mxArrayExtractStringField(matlabArray, "type");
  if(mxType!=type)
    cerr << "Error mismatch between saved type, " << mxType << ", and Class type, " << type << "." << endl;
  mxArray* transformArray = mxArrayExtractMxArrayField(matlabArray, "transforms");
  // transforms field.
  transformsFromMxArray(transformArray);
    
  // TODO priors ... need to deal with priors
  extractParamFromMxArray(matlabArray);
}
void CKern::getGradPrior(CMatrix& G) const
{
}
void CKern::getPriorLogProb(CMatrix& G) const
{
}
void CKern::extractParamFromMxArray(const mxArray* matlabArray)
{
  inputDim = mxArrayExtractIntField(matlabArray, "inputDimension");
  nParams = mxArrayExtractIntField(matlabArray, "nParams");
  string pName;
  for(int i=0; i<nParams; i++)
    {
      pName=getParamName(i);
      setParam(mxArrayExtractDoubleField(matlabArray, pName), i);
    }
}
void CKern::addParamToMxArray(mxArray* matlabArray) const
{
  mxAddField(matlabArray, "nParams");
  mxSetField(matlabArray, 0, "nParams", convertMxArray((double)nParams));
  string pName;
  for(int i=0; i<nParams; i++)
    {
      pName = getParamName(i);
      mxAddField(matlabArray, pName.c_str());      
      mxSetField(matlabArray, 0, pName.c_str(), convertMxArray(getParam(i))); 
    }
}

void CKern::getGradTransParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd) const
{
  assert(g.getRows()==1);
  assert(g.getCols()==getNumParams());
  getGradParams(g, X, cvGrd);
  double val;
  double param;
  for(int i=0; i<getNumTransforms(); i++)
    {
      val=g.getVal(getTransformIndex(i));
      param=getParam(getTransformIndex(i));
      g.setVal(val*getTransformGradFact(param, i), getTransformIndex(i));
    }  
}
bool CKern::equals(const CKern& kern, const double tol) const
{
  if(getType()!=kern.getType())
    return false;
  if(getNumParams()!=kern.getNumParams())
    return false;
  CMatrix params(1, getNumParams());
  getParams(params);
  CMatrix kernParams(1, getNumParams());
  kern.getParams(kernParams);
  if(!params.equals(kernParams, tol))
    return false;
  return true;
}

// The Compound kernel.
CCmpndKern::CCmpndKern()
{
}
CCmpndKern::CCmpndKern(const int inDim)
{
  initialiseKern(inDim);
  setInitParam();
}
CCmpndKern::CCmpndKern(const CMatrix& X) 
{
  initialiseKern(X);
  setInitParam();
}  
// create a Ckern given a matlab kern structure.
CCmpndKern::CCmpndKern(mxArray* kern)
{
}
// Class destructor
CCmpndKern::~CCmpndKern()
{
}
CCmpndKern::CCmpndKern(const CCmpndKern& kern) : components(kern.components)
{
  initialiseKern(kern);
  setInitParam();
  for(int i=0; i<components.size(); i++)
    addKern(components[i]->clone()); 
  //  components(kern.components);// = kern.components;
}
void CCmpndKern::setInitParam()
{
  setType("cmpnd");
  setName("compound");
  nParams=0;
}
double CCmpndKern::diagComputeElement(const CMatrix& X, const int index) const
{
  double y=0.0;
  for(int i=0; i<components.size(); i++)
    y+=components[i]->diagComputeElement(X, index);
  return y;
}
void CCmpndKern::diagCompute(CMatrix& d, const CMatrix& X) const
{
  assert(d.getCols()==1);
  assert(X.rowsMatch(d));
  d.zeros();
  CMatrix dStore(d.getRows(), d.getCols(), 0.0);
  for(int i=0; i < components.size(); i++)
    {
      components[i]->diagCompute(dStore, X);
      d.axpy(dStore, 1.0);
    }
}
void CCmpndKern::setParam(const double val, const int paramNo)
{
  int start = 0;
  int end = 0;
  for(int i=0; i<components.size(); i++)
    {
      end = start+components[i]->getNumParams()-1;
      if(paramNo <= end)
	{
	  components[i]->setParam(val, paramNo-start);
	  return;
	}
      
      start = end + 1;
    }
}
// Parameters are kernel parameters
double CCmpndKern::getParam(int paramNo) const
{
  int start = 0;
  int end = 0;
  for(int i=0; i<components.size(); i++)
    {
      end = start+components[i]->getNumParams()-1;
      if(paramNo <= end)
	return components[i]->getParam(paramNo-start);
      start = end + 1;
    }
}
string CCmpndKern::getParamName(int paramNo) const
{
  int start = 0;
  int end = 0;
  for(int i=0; i<components.size(); i++)
    {
      end = start+components[i]->getNumParams()-1;
      if(paramNo <= end)
	return components[i]->getType() + components[i]->getParamName(paramNo-start);
      start = end + 1;
    }
}
void CCmpndKern::getGradX(vector<CMatrix*> G, const CMatrix& X, const CMatrix& X2) const
{
  vector<CMatrix>::const_iterator iter2;
//   //vector<CMatrix> G;
//   //CMatrix g(X2.getRows(), X2.getCols(), 0.0);
//   for(int i=0; i<components.size(); i++)
//     {
//       vector<CMatrix> subG = components[i]->getGradX(X, X2);
//       for(iter2=subG.begin(); iter2!=subG.end(); ++iter2)
// 	{
// 	  g+=*iter2;
// 	}	
//       G.push_back(g);
//     }
//   return G;
}
void CCmpndKern::getDiagGradX(CMatrix& G, const CMatrix& X) const
{
  assert(G.dimensionsMatch(X));
  G.zeros();
  CMatrix tempG(G.getRows(), G.getCols());
  for(int i=0; i<components.size(); i++)
    {
      components[i]->getDiagGradX(tempG, X);
      G.axpy(tempG, 1.0);
    }
}
double CCmpndKern::getDiagGradXElement(const CMatrix& X, const int i, const int j) const
{
  double g = 0.0;
  for(int comp=0; i<components.size(); comp++)
    {
      g+=components[comp]->getDiagGradXElement(X, i, j);
    }
  return g;
}  
double CCmpndKern::getWhite() const
{
  double white = 0.0;
  for(int i=0; i<components.size(); i++)
    white += components[i]->getWhite();
  return white;
}

double CCmpndKern::computeElement(const CMatrix& X1, const int index1, 
			   const CMatrix& X2, const int index2) const
{
  double y=0.0;
  for(int i=0; i<components.size(); i++)
    y+=components[i]->computeElement(X1, index1, X2, index2);
  return y;
}

void CCmpndKern::compute(CMatrix& K, const CMatrix& X) const
{
  assert(K.rowsMatch(X));
  assert(K.isSquare());
  K.zeros();
  CMatrix K2(K.getRows(), K.getCols());
  for(int i=0; i<components.size(); i++)
    {
      components[i]->compute(K2, X);
      K.axpy(K2, 1.0);
    }
}

void CCmpndKern::compute(CMatrix& K, const CMatrix& X, const CMatrix& X2) const
{
  assert(K.rowsMatch(X));
  assert(K.getCols()==X2.getRows());
  K.zeros();
  CMatrix K2(K.getRows(), K.getCols());
  for(int i=0; i<components.size(); i++)
    {
      components[i]->compute(K2, X, X2);
      K.axpy(K2, 1.0);
    }
}
void CCmpndKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& covGrad) const
{
  assert(g.getRows()==1);
  assert(g.getCols()==getNumParams());
  int start = 0;
  int end = 0;
  for(int i=0; i<components.size(); i++)
    {
      end = start+components[i]->getNumParams()-1;
      CMatrix subg(1, components[i]->getNumParams());
      components[i]->getGradParams(subg, X, covGrad);
      g.setMatrix(0, start, subg);
      start = end+1;
    }
  addPriorGrad(g);
}
double CCmpndKern::getGradParam(const int index, const CMatrix& X, const CMatrix& covGrad) const
{
  assert(index>=0);
  assert(index<nParams);
  int start=0;
  int end=0;
  for(int i=0; i<components.size(); i++)
    {
      end = start+components[i]->getNumParams()-1;
      if(index<end)
	return components[i]->getGradParam(index-start, X, covGrad);
    }
}
void CCmpndKern::extractParamFromMxArray(const mxArray* matlabArray)
{
  nParams=0;
  mxArray* compArray=mxArrayExtractMxArrayField(matlabArray, "comp");
  int numComps = mxGetNumberOfElements(compArray);
  mxArray* compElement;
  string kernType;
  int kernNum = 0;
  for(int i=0; i<numComps; i++)
    {
      // TODO this is a pretty ugly way of doing things ... perhaps dynamic casting can be used?
      compElement = mxGetCell(compArray, i);
      kernType = mxArrayExtractStringField(compElement, "type");
      if(kernType == "lin")
	kernNum = addKern(new CLinKern(inputDim));
      else if(kernType == "rbf")  
	kernNum = addKern(new CRbfKern(inputDim));
      else if(kernType == "white")
	kernNum = addKern(new CWhiteKern(inputDim));
      else if(kernType == "bias")
	kernNum = addKern(new CBiasKern(inputDim));
      else if(kernType == "cmpnd")
	kernNum = addKern(new CCmpndKern(inputDim));
      
      components[kernNum]->fromMxArray(compElement);
    }
}
int CCmpndKern::addKern(CKern* kern)
{
  components.push_back(kern);
  int oldNParams = nParams;
  nParams+=kern->getNumParams();
  for(int i=0; i<kern->getNumTransforms(); i++)
    addTransform(kern->getTransform(i), kern->getTransformIndex(i)+oldNParams);      
  return components.size()-1;
}
void CCmpndKern::addParamToMxArray(mxArray* matlabArray) const
{
  // Add comp field to mxArray.
  int dims[0];
  mxAddField(matlabArray, "comp");
  dims[0] = components.size();
  mxArray* compArray = mxCreateCellArray(1, dims);
  mxSetField(matlabArray, 0, "comp", 
	     compArray); 
  mxArray* currentKern;
  for(int i=0; i<components.size(); i++)
    {
      currentKern = components[i]->toMxArray();
      mxAddField(currentKern, "index");
      mxSetCell(compArray, i, currentKern);
    }

}
// the white noise kernel.
CWhiteKern::CWhiteKern()
{
}
CWhiteKern::CWhiteKern(const int inDim)
{
  initialiseKern(inDim);
  setInitParam();
}
CWhiteKern::CWhiteKern(const CMatrix& X)
{
  initialiseKern(X);
  setInitParam();
}  
// create a Ckern given a matlab kern structure.
CWhiteKern::CWhiteKern(mxArray* kern)
{
}
  // Class destructor
CWhiteKern::~CWhiteKern()
{
}
CWhiteKern::CWhiteKern(const CWhiteKern& kern)
{
  initialiseKern(kern);
  setInitParam();
  variance = kern.variance;
}
void CWhiteKern::setInitParam()
{
  nParams = 1;
  setType("white");
  setName("white noise");
  setParamName("variance", 0);
  variance = exp(-2.0);
  addTransform(new CNegLogLogitTransform, 0);
  // TODO Add prior functionality priors.push_back(new CDist);
}

inline double CWhiteKern::diagComputeElement(const CMatrix& X, const int index) const
{
  return variance;
}
void CWhiteKern::diagCompute(CMatrix& d, const CMatrix& X) const
{
  assert(d.getCols()==1);
  assert(X.rowsMatch(d));
  d.setVals(variance);
}
void CWhiteKern::setParam(const double val, const int paramNo)
{
  assert(paramNo==0);
  switch(paramNo)
    {
    case 0:
      variance = val;
      break;
    default:
      cerr<<"Parameter doesn't exist.";
    }
}
// Parameters are kernel parameters
double CWhiteKern::getParam(const int paramNo) const
{
  assert(paramNo==0);
  switch(paramNo)
    {
    case 0:
      return variance;
      break;
    default:
      cerr<<"Parameter doesn't exist.";
    }
}
void CWhiteKern::getGradX(vector<CMatrix*> G, const CMatrix& X, const CMatrix& X2) const
{
  assert(G.size()==X.getRows());
  for(int i=0; i<G.size(); i++)
    {
      assert(G[i]->getRows()==X.getRows());
      assert(G[i]->getCols()==X.getCols());
      G[i]->zeros();
    }
}
void CWhiteKern::getDiagGradX(CMatrix& G, const CMatrix& X) const
{
  assert(G.dimensionsMatch(X));
  G.zeros();
}
inline double CWhiteKern::getDiagGradXElement(const CMatrix& X, const int i, const int j) const
{
  return 0.0;
}
double CWhiteKern::getWhite() const
{
  return variance;
}

inline double CWhiteKern::computeElement(const CMatrix& X1, const int index1,
				  const CMatrix& X2, const int index2) const
{
  return 0.0;
}

void CWhiteKern::compute(CMatrix& K, const CMatrix& X) const
{
  assert(K.rowsMatch(X));
  assert(K.isSquare());
  K.zeros();
  for(int i=0; i<K.getRows(); i++)
    K.setVal(variance, i, i);
  K.setSymmetric(true);
}

void CWhiteKern::compute(CMatrix& K, const CMatrix& X, const CMatrix& X2) const
{
  assert(K.rowsMatch(X));
  assert(K.getCols()==X2.getRows());
  K.zeros();
}
double CWhiteKern::getGradParam(const int index, const CMatrix& X, const CMatrix& covGrad) const
{
  assert(index==0);
  return trace(covGrad);
}

// the bias kernel.
CBiasKern::CBiasKern()
{
}
CBiasKern::CBiasKern(const int inDim)
{
  initialiseKern(inDim);
  setInitParam();
}
CBiasKern::CBiasKern(const CMatrix& X)
{
  initialiseKern(X);
  setInitParam();
}  
// create a Ckern given a matlab kern structure.
CBiasKern::CBiasKern(mxArray* kern)
{
}
  // Class destructor
CBiasKern::~CBiasKern()
{
}
CBiasKern::CBiasKern(const CBiasKern& kern)
{
  initialiseKern(kern);
  setInitParam();
  variance = kern.variance;
  
}
void CBiasKern::setInitParam()
{
  nParams = 1;
  setType("bias");
  setName("bias");
  setParamName("variance", 0);
  variance = exp(-2.0);
  addTransform(new CNegLogLogitTransform, 0);
  // TODO Add prior functionality priors.push_back(new CDist);
}

double CBiasKern::diagComputeElement(const CMatrix& X, const int index) const
{
  return variance;
}
void CBiasKern::diagCompute(CMatrix& d, const CMatrix& X) const
{
  assert(d.getCols()==1);
  assert(X.rowsMatch(d));
  d.setVals(variance);
}
// Parameters are kernel parameters
void CBiasKern::setParam(const double val, const int paramNo)
{
  assert(paramNo==0);
  switch(paramNo)
    {
    case 0:
      variance = val;
      break;
    default:
      cerr<<"Parameter doesn't exist.";
    }
}
double CBiasKern::getParam(const int paramNo) const
{
  assert(paramNo==0);
  switch(paramNo)
    {
    case 0:
      return variance;
      break;
    default:
      cerr<<"Parameter doesn't exist.";
    }
}
void CBiasKern::getGradX(vector<CMatrix*> G, const CMatrix& X, const CMatrix& X2) const
{
  assert(G.size()==X.getRows());
  for(int i=0; i<G.size(); i++)
    {
      assert(G[i]->getRows()==X.getRows());
      assert(G[i]->getCols()==X.getCols());
      G[i]->zeros();
    }
}
void CBiasKern::getDiagGradX(CMatrix& G, const CMatrix& X) const
{
  assert(G.dimensionsMatch(X));
  G.zeros();
}
inline double CBiasKern::getDiagGradXElement(const CMatrix& X, const int i, const int j) const
{
  return 0.0;
}
double CBiasKern::getWhite() const
{
  return 0.0;
}

inline double CBiasKern::computeElement(const CMatrix& X1, const int index1, 
				 const CMatrix& X2, const int index2) const
{
  return variance;
}

void CBiasKern::compute(CMatrix& K, const CMatrix& X) const
{
  assert(K.rowsMatch(X));
  assert(K.isSquare());
  K.setVals(variance);
  K.setSymmetric(true);
}

void CBiasKern::compute(CMatrix& K, const CMatrix& X, const CMatrix& X2) const
{
  assert(K.rowsMatch(X));
  assert(K.getCols()==X2.getRows());
  K.setVals(variance);
}
double CBiasKern::getGradParam(const int index, const CMatrix& X, const CMatrix& covGrad) const 
{
  assert(index==0);
  return sum(covGrad);
}
// the RBF kernel.
CRbfKern::CRbfKern()
{
}
CRbfKern::CRbfKern(const int inDim)
{
  initialiseKern(inDim);
  setInitParam();
}
CRbfKern::CRbfKern(const CMatrix& X)
{
  initialiseKern(X);
  setInitParam();
}  
// create a Ckern given a matlab kern structure.
CRbfKern::CRbfKern(mxArray* kern)
{
}
  // Class destructor
CRbfKern::~CRbfKern()
{
}
CRbfKern::CRbfKern(const CRbfKern& kern)
{
  initialiseKern(kern);
  setInitParam();
  variance = kern.variance;
  inverseWidth = inverseWidth;
}
void CRbfKern::setInitParam()
{
  nParams = 2;
  setType("rbf");
  setName("RBF");
  setParamName("inverseWidth", 0);
  inverseWidth = 1.0;
  setParamName("variance", 1);
  variance = 1.0;
  addTransform(new CNegLogLogitTransform, 0);
  addTransform(new CNegLogLogitTransform, 1);
  // TODO Add prior functionality priors.push_back(new CDist);
}

inline double CRbfKern::diagComputeElement(const CMatrix& X, const int index) const
{
  return variance;
}
void CRbfKern::diagCompute(CMatrix& d, const CMatrix& X) const
{
  assert(d.getCols()==1);
  assert(X.rowsMatch(d));
  d.setVals(variance);
}
// Parameters are kernel parameters
void CRbfKern::setParam(const double val, const int paramNo)
{
  assert(paramNo < nParams);
  switch(paramNo)
    {
    case 0:
      inverseWidth = val;
      break;
    case 1:
      variance = val;
      break;
    default:
      cerr << "Parameter doesn't exist.";
    }
}
double CRbfKern::getParam(const int paramNo) const
{
  assert(paramNo < nParams);
  switch(paramNo)
    {
    case 0:
      return inverseWidth;
      break;
    case 1:
      return variance;
      break;
    default:
      cerr << "Parameter doesn't exist.";
    }
}

void CRbfKern::getGradX(vector<CMatrix*> G, const CMatrix& X, const CMatrix& X2) const
{
  // TODO this needs fixing.
  assert(G.size()==X.getRows());
  for(int i=0; i<G.size(); i++)
    {
      assert(G[i]->getRows()==X.getRows());
      assert(G[i]->getCols()==X.getCols());
      G[i]->zeros();
    }
}
void CRbfKern::getDiagGradX(CMatrix& G, const CMatrix& X) const
{
  assert(G.dimensionsMatch(X));
  G.zeros();
}
double CRbfKern::getDiagGradXElement(const CMatrix& X, const int i, const int j) const
{
  return 0.0;
}
double CRbfKern::getWhite() const
{
  return 0.0;
}

double CRbfKern::computeElement(const CMatrix& X1, const int index1, 
			 const CMatrix& X2, const int index2) const
{
  double k = X1.dist2Row(index1, X2, index2);
  k = 0.5*k*inverseWidth;
  k = variance*exp(-k);
  return k;
}
void CRbfKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& covGrad) const
{
  assert(g.getRows()==1);
  assert(g.getCols()==nParams);
  double dist2;
  double g1=0.0;
  double g2=0.0;
  double k;
  for(int i=0; i<X.getRows(); i++)
    for(int j=0; j<X.getRows(); j++)
      {
	dist2 = X.dist2Row(i, X, j);
	k = exp(-dist2*0.5*inverseWidth);
	g1 -= 0.5*k*variance*dist2*covGrad.getVal(i, j);
	g2 += k*covGrad.getVal(i, j);
      }
  g.setVal(g1, 0);
  g.setVal(g2, 1);
  addPriorGrad(g);
}
double CRbfKern::getGradParam(const int index, const CMatrix& X, const CMatrix& covGrad) const
{
  assert(index>=0);
  assert(index<nParams);
//   switch(index)
//     {
//     case 0:
//       // todo 

//       break;
//     case 1:
//       // todo
//       break;
//     case default:
//       assert(0);
//     }

//   //  CMatrix KcovGrad = 
//   // CMatrix K = compute(X);
//   return sum(covGrad);
}


// the Linear kernel.
CLinKern::CLinKern()
{
}
CLinKern::CLinKern(const int inDim)
{
  initialiseKern(inDim);
  setInitParam();
}
CLinKern::CLinKern(const CMatrix& X)
{
  initialiseKern(X);
  setInitParam();
}  
// create a Ckern given a matlab kern structure.
CLinKern::CLinKern(mxArray* kern)
{
}
  // Class destructor
CLinKern::~CLinKern()
{
}
CLinKern::CLinKern(const CLinKern& kern)
{
  initialiseKern(kern);
  setInitParam();
  variance = kern.variance;
}
void CLinKern::setInitParam()
{
  nParams = 1;
  setType("lin");
  setName("linear");
  setParamName("variance", 0);
  variance = 1.0;
  addTransform(new CNegLogLogitTransform, 0);
  // TODO add prior functionality priors.push_back(new CDist);
}

double CLinKern::diagComputeElement(const CMatrix& X, const int index1) const
{
  return variance*X.norm2Row(index1);  
}
// Parameters are kernel parameters
void CLinKern::setParam(const double val, const int paramNo)
{
  assert(paramNo==0);
  switch(paramNo)
    {
    case 0:
      variance = val;
      break;
    default:
      cerr << "Parameter doesn't exist.";
    }
}
double CLinKern::getParam(const int paramNo) const
{
  assert(paramNo==0);
  switch(paramNo)
    {
    case 0:
      return variance;
      break;
    default:
      cerr << "Parameter doesn't exist.";
    }
}
void CLinKern::getGradX(vector<CMatrix*> G, const CMatrix& X, const CMatrix& X2) const
{
  // TODO This needs fixing.
  assert(G.size()==X.getRows());
  for(int i=0; i<G.size(); i++)
    {
      assert(G[i]->getRows()==X.getRows());
      assert(G[i]->getCols()==X.getCols());
      G[i]->zeros();
    }
}
void CLinKern::getDiagGradX(CMatrix& G, const CMatrix& X) const
{
  assert(G.dimensionsMatch(X));
  G.deepCopy(X);
  G.scale(2.0*variance);
}
double CLinKern::getDiagGradXElement(const CMatrix& X, const int i, const int j) const
{
  double g = X.getVal(i, j);
  g*=2.0*variance;
  return g;
}
double CLinKern::getWhite() const
{
  return 0.0;
}

double CLinKern::computeElement(const CMatrix& X1, const int index1, 
			  const CMatrix& X2, const int index2) const
{
  return variance*X1.dotRowRow(index1, X2, index2);
}

void CLinKern::compute(CMatrix& K, const CMatrix& X) const
{
  assert(K.rowsMatch(X));
  assert(K.isSquare());
  K.setSymmetric(true);
  K.syrk(X, variance, 0.0, "u", "n");
}

void CLinKern::compute(CMatrix& K, const CMatrix& X, const CMatrix& X2) const
{
  assert(K.rowsMatch(X));
  assert(K.getCols()==X2.getRows());
  K.gemm(X, X2, variance, 0.0, "n", "t");
}
double CLinKern::getGradParam(const int index, const CMatrix& X, const CMatrix& covGrad) const
{
  assert(index==0);
  assert(X.rowsMatch(covGrad));
  assert(covGrad.isSquare());
  double dot;
  double g1=0.0;
  for(int i=0; i<X.getRows(); i++)
    for(int j=0; j<X.getRows(); j++)
      {
	dot = X.dotRowRow(i, X, j);
	g1 += dot*covGrad.getVal(i, j);
      }
  return g1;
}
ostream& operator<<(ostream& out, const CKern& kern)
{
  out <<  kern.display(out);
  return out;
}
