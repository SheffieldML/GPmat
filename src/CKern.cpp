#include "CKern.h"

using namespace std;

ostream& CKern::display(ostream& os) const
{
  os << getName() << " kernel:" << endl;
  for(unsigned int i=0; i<nParams; i++)
  {
    os << getParamName(i) << ": " << getParam(i) << endl;
  }
  return os;
}

void CKern::writeParamsToStream(ostream& out) const
{
  writeToStream(out, "baseType", getBaseType());
  writeToStream(out, "type", getType());
  writeToStream(out, "inputDim", getInputDim());
  writeToStream(out, "numParams", getNumParams());
  CMatrix par(1, getNumParams());
  getParams(par);
  par.toStream(out);
  writeToStream(out, "numPriors", getNumPriors());
  writePriorsToStream(out);
}

void CKern::getGradPrior(CMatrix& G) const
{
}

void CKern::getPriorLogProb(CMatrix& G) const
{
}

void CKern::getGradTransParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==getNumParams());
  getGradParams(g, X, X2, cvGrd, regularise);
  double val;
  double param;
  for(unsigned int i=0; i<getNumTransforms(); i++)
  {
    val=g.getVal(getTransformIndex(i));
    param=getParam(getTransformIndex(i));
    g.setVal(val*getTransformGradFact(param, i), getTransformIndex(i));
  }
}
void CKern::getGradTransParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==getNumParams());
  getGradParams(g, X, cvGrd, regularise);
  double val;
  double param;
  for(unsigned int i=0; i<getNumTransforms(); i++)
  {
    val=g.getVal(getTransformIndex(i));
    param=getParam(getTransformIndex(i));
    g.setVal(val*getTransformGradFact(param, i), getTransformIndex(i));
  }  
}
void CKern::getDiagGradTransParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==getNumParams());
  getDiagGradParams(g, X, cvGrd, regularise);
  double val;
  double param;
  for(unsigned int i=0; i<getNumTransforms(); i++)
  {
    val=g.getVal(getTransformIndex(i));
    param=getParam(getTransformIndex(i));
    g.setVal(val*getTransformGradFact(param, i), getTransformIndex(i));
  }  
}
bool CKern::equals(const CKern& kern, double tol) const
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
// The Component kernel
void CComponentKern::readParamsFromStream(istream& in) 
{
//   string tbaseType = getBaseTypeStream(in);
//   if(tbaseType != getBaseType())
//     throw ndlexceptions::StreamFormatError("baseType", "Error mismatch between saved base type, " + tbaseType + ", and Class base type, " + getType() + ".");
//   string ttype = getTypeStream(in);
//   if(ttype != getType())
//     throw ndlexceptions::StreamFormatError("type", "Error mismatch between saved type, " + ttype + ", and Class type, " + getType() + ".");
  setInputDim(readIntFromStream(in, "inputDim"));
  unsigned int nPar = readIntFromStream(in, "numParams");
  unsigned int numKerns = readIntFromStream(in, "numKerns");
  for(unsigned int i=0; i<numKerns; i++)
  {
    if(i<components.size())
      components[i]->fromStream(in);
    else
      addKern(readKernFromStream(in));
  }
  
}
void CComponentKern::writeParamsToStream(ostream& out) const
{
  writeToStream(out, "baseType", getBaseType());
  writeToStream(out, "type", getType());
  writeToStream(out, "inputDim", getInputDim());
  writeToStream(out, "numParams", getNumParams());
  writeToStream(out, "numKerns", (unsigned int)components.size());
  for(unsigned int i=0; i<components.size(); i++)
  {
    components[i]->toStream(out);
  }
}


// The Compound kernel.
CCmpndKern::CCmpndKern() : CComponentKern()
{
  _init();
}
CCmpndKern::CCmpndKern(unsigned int inDim) : CComponentKern(inDim)
{
  _init();
  setInputDim(inDim);
}
CCmpndKern::CCmpndKern(const CMatrix& X) : CComponentKern(X)
{
  _init();
  setInputDim(X.getCols());
}  
CCmpndKern::CCmpndKern(const CCmpndKern& kern) : CComponentKern(kern)
{
  _init();
  setInputDim(kern.getInputDim());
  for(size_t i=0; i<components.size(); i++)
    addKern(components[i]->clone()); 
}
// Class destructor
CCmpndKern::~CCmpndKern()
{
  for(size_t i=0; i<components.size(); i++)
    delete components[i];
}
void CCmpndKern::_init()
{
  nParams=0;
  setType("cmpnd");
  setName("compound");
  setStationary(true);
}
void CCmpndKern::setInitParam()
{
}
double CCmpndKern::diagComputeElement(const CMatrix& X, unsigned int index) const
{
  double y=0.0;
  for(size_t i=0; i<components.size(); i++)
    y+=components[i]->diagComputeElement(X, index);
  return y;
}
void CCmpndKern::diagCompute(CMatrix& d, const CMatrix& X) const
{
  DIMENSIONMATCH(d.getCols()==1);
  DIMENSIONMATCH(X.rowsMatch(d));
  d.zeros();
  CMatrix dStore(d.getRows(), d.getCols(), 0.0);
  for(size_t i=0; i < components.size(); i++)
  {
    components[i]->diagCompute(dStore, X);
    d.axpy(dStore, 1.0);
  }
}
void CCmpndKern::getGradX(CMatrix& gX, const CMatrix& X, unsigned int row, const CMatrix& X2, bool addG) const
{
  DIMENSIONMATCH(gX.getRows() == X2.getRows());
  BOUNDCHECK(row < X.getRows());
  DIMENSIONMATCH(X.getCols()==X2.getCols());
  if(!addG)
    gX.zeros();
  for(size_t i=0; i<components.size(); i++)
    components[i]->getGradX(gX, X, row, X2, true);
}
void CCmpndKern::getDiagGradX(CMatrix& gX, const CMatrix& X, bool addG) const
{
  DIMENSIONMATCH(gX.dimensionsMatch(X));
  if(!addG)
    gX.zeros();
  for(size_t i=0; i<components.size(); i++)
  {
    components[i]->getDiagGradX(gX, X, true);
  }
}
double CCmpndKern::getWhite() const
{
  double white = 0.0;
  for(size_t i=0; i<components.size(); i++)
    white += components[i]->getWhite();
  return white;
}
double CCmpndKern::getVariance() const
{
  double totVariance = 0.0;
  for(size_t i=0; i<components.size(); i++)
    totVariance += components[i]->getVariance();
  return totVariance;
}

double CCmpndKern::computeElement(const CMatrix& X1, unsigned int index1, 
				  const CMatrix& X2, unsigned int index2) const
{
  double y=0.0;
  for(size_t i=0; i<components.size(); i++)
    y+=components[i]->computeElement(X1, index1, X2, index2);
  return y;
}

// void CCmpndKern::compute(CMatrix& K, const CMatrix& X) const
// {
//   DIMENSIONMATCH(K.rowsMatch(X));
//   MATRIXPROPERTIES(K.isSquare());
//   K.zeros();
//   CMatrix K2(K.getRows(), K.getCols());
//   for(size_t i=0; i<components.size(); i++)
//   {
//     components[i]->compute(K2, X);
//     K.axpy(K2, 1.0);
//   }
// }
// void CCmpndKern::compute(CMatrix& K, const CMatrix& X, const CMatrix& X2, unsigned int row) const
// {
//   DIMENSIONMATCH(K.rowsMatch(X));
//   DIMENSIONMATCH(K.getCols()==1);
//   //CMatrix K2(K.getRows(), K.getCols());
//   for(unsigned int i=0; i<X.getRows(); i++)
//   {
//     double kval=0.0;
//     for(unsigned int k=0; k<components.size(); k++)
//       kval+=components[k]->computeElement(X, i, X2, row);
//     K.setVal(kval, i, row);
//   }
// }
// void CCmpndKern::compute(CMatrix& K, const CMatrix& X, const CMatrix& X2) const
// {
//   DIMENSIONMATCH(K.rowsMatch(X));
//   DIMENSIONMATCH(K.getCols()==X2.getRows());
//   //CMatrix K2(K.getRows(), K.getCols());
//   for(unsigned int i=0; i<X.getRows(); i++)
//   {
//     for(unsigned int j=0; j<X2.getRows(); j++)
//     {
//       double kval=0.0;
//       for(unsigned int k=0; k<components.size(); k++)
// 	kval+=components[k]->computeElement(X, i, X2, j);
//       K.setVal(kval, i, j);
//     }
//   }
// }
void CCmpndKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==getNumParams());
  unsigned int start = 0;
  unsigned int end = 0;
  for(size_t i=0; i<components.size(); i++)
  {
    end = start+components[i]->getNumParams()-1;
    CMatrix subg(1, components[i]->getNumParams());
    components[i]->getGradParams(subg, X, X2, covGrad, regularise);
    g.setMatrix(0, start, subg);
    start = end+1;
  }
}
void CCmpndKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& covGrad, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==getNumParams());
  unsigned int start = 0;
  unsigned int end = 0;
  for(size_t i=0; i<components.size(); i++)
  {
    end = start+components[i]->getNumParams()-1;
    CMatrix subg(1, components[i]->getNumParams());
    components[i]->getGradParams(subg, X, covGrad, regularise);
    g.setMatrix(0, start, subg);
    start = end+1;
  }
}
double CCmpndKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad) const
{
  
  BOUNDCHECK(index<nParams);
  unsigned int start=0;
  unsigned int end=0;
  for(size_t i=0; i<components.size(); i++)
  {
    end = start+components[i]->getNumParams()-1;
    if(index<end)
      return components[i]->getGradParam(index-start, X, X2, covGrad);
    start = end + 1;
  }
  return -1;
}
double CCmpndKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& covGrad) const
{
  
  BOUNDCHECK(index<nParams);
  unsigned int start=0;
  unsigned int end=0;
  for(size_t i=0; i<components.size(); i++)
  {
    end = start+components[i]->getNumParams()-1;
    if(index<end)
      return components[i]->getGradParam(index-start, X, covGrad);
    start = end + 1;
  }
  return -1;
}

// The Tensor kernel.
CTensorKern::CTensorKern() : CComponentKern()
{
  _init();
}
CTensorKern::CTensorKern(unsigned int inDim) : CComponentKern(inDim)
{
  _init();
  setInputDim(inDim);
}
CTensorKern::CTensorKern(const CMatrix& X) : CComponentKern(X)
{
  _init();
  setInputDim(X.getCols());
}  
CTensorKern::CTensorKern(const CTensorKern& kern) : CComponentKern(kern)
{

  _init();
  setInputDim(kern.getInputDim());
  for(size_t i=0; i<components.size(); i++)
  {
    addKern(components[i]->clone()); 
  }
}
CTensorKern::CTensorKern(const CTensorKern& kern, unsigned int comp) : CComponentKern(kern)
{
  _init();
  BOUNDCHECK(comp<kern.components.size());
  setInputDim(kern.getInputDim());
  setInitParam();
  nParams = kern.getNumParams();
  nParams -= kern.components[comp]->getNumParams();
  for(size_t i=0; i<comp; i++)
  {
    addKern(kern.components[i]->clone());
  }
  for(size_t i=comp+1; i<kern.components.size(); i++)
  {
    addKern(kern.components[i]->clone());
  }
}
// Class destructor
CTensorKern::~CTensorKern()
{
  for(size_t i=0; i<components.size(); i++)
    delete components[i];

}
void CTensorKern::_init()
{
  nParams=0;
  setType("tensor");
  setName("tensor");
  setStationary(true);
}
void CTensorKern::setInitParam()
{
}
unsigned int CTensorKern::addKern(const CKern* kern)
{
  if(kern->getWhite())
  {
    throw ndlexceptions::Error("Components of a tensor kernel should not contain white noise");

  }
  return CComponentKern::addKern(kern);
}
double CTensorKern::diagComputeElement(const CMatrix& X, unsigned int index) const
{
  double y=1.0;
  for(size_t i=0; i<components.size(); i++)
    y*=components[i]->diagComputeElement(X, index);
  return y;
}
void CTensorKern::diagCompute(CMatrix& d, const CMatrix& X) const
{
  DIMENSIONMATCH(d.getCols()==1);
  DIMENSIONMATCH(X.rowsMatch(d));
  d.ones();
  CMatrix dStore(d.getRows(), d.getCols(), 0.0);
  for(size_t i=0; i < components.size(); i++)
  {
    components[i]->diagCompute(dStore, X);
    d.dotMultiply(dStore);
  }
}
double CTensorKern::getWhite() const
{
  return 0.0;
}
double CTensorKern::getVariance() const
{
  double totVariance = 1.0;
  for(size_t i=0; i<components.size(); i++)
    totVariance *= components[i]->getVariance();
  return totVariance;
}

double CTensorKern::computeElement(const CMatrix& X1, unsigned int index1, 
				   const CMatrix& X2, unsigned int index2) const
{
  double y=1.0;
  for(size_t i=0; i<components.size(); i++)
    y*=components[i]->computeElement(X1, index1, X2, index2);
  return y;
}

// void CTensorKern::compute(CMatrix& K, const CMatrix& X) const
// {
//   DIMENSIONMATCH(K.rowsMatch(X));
//   DIMENSIONMATCH(K.isSquare());
//   K.ones();
//   CMatrix K2(K.getRows(), K.getCols());
//   for(size_t i=0; i<components.size(); i++)
//   {
//     components[i]->compute(K2, X);
//     K.dotMultiply(K2);
//   }
// }

// void CTensorKern::compute(CMatrix& K, const CMatrix& X, const CMatrix& X2) const
// {
//   DIMENSIONMATCH(K.rowsMatch(X));
//   DIMENSIONMATCH(K.getCols()==X2.getRows());
//   //CMatrix K2(K.getRows(), K.getCols());
//   for(unsigned int i=0; i<X.getRows(); i++)
//   {
//     for(unsigned int j=0; j<X2.getRows(); j++)
//     {
//       double kval=1.0;
//       for(size_t k=0; k<components.size(); k++)
// 	kval*=components[k]->computeElement(X, i, X2, j);
//       K.setVal(kval, i, j);
//     }
//   }
// }
// void CTensorKern::compute(CMatrix& K, const CMatrix& X, const CMatrix& X2, unsigned int row) const
// {
//   DIMENSIONMATCH(K.rowsMatch(X));
//   DIMENSIONMATCH(K.getCols()==1);
//   //CMatrix K2(K.getRows(), K.getCols());
//   for(unsigned int i=0; i<X.getRows(); i++)
//   {
//     double kval=1.0;
//     for(size_t k=0; k<components.size(); k++)
//       kval*=components[k]->computeElement(X, i, X2, row);
//     K.setVal(kval, i, row);
//   }
// }
void CTensorKern::getGradX(CMatrix& gX, const CMatrix& X, unsigned int row, const CMatrix& X2, bool addG) const
{
  DIMENSIONMATCH(gX.getRows() == X2.getRows());
  BOUNDCHECK(row < X.getRows());
  DIMENSIONMATCH(X.getCols()==X2.getCols());
  if(!addG)
    gX.zeros();
  CMatrix gXTemp(gX.getRows(), gX.getCols());
  CMatrix Kslash(X2.getRows(), 1);
  for(size_t i=0; i<components.size(); i++)
  {
    CTensorKern* slashKern = new CTensorKern(*this, i);
    components[i]->getGradX(gXTemp, X, row, X2);
    slashKern->compute(Kslash, X2, X, row);
    for(unsigned int j=0; j<getInputDim(); j++)
    {
      gXTemp.dotMultiplyColCol(j, Kslash, 0);
    }
    gX.add(gXTemp);
  }
}
void CTensorKern::getDiagGradX(CMatrix& gX, const CMatrix& X, bool addG) const
{
  DIMENSIONMATCH(gX.dimensionsMatch(X));
  if(!addG)
    gX.zeros();
  CMatrix gXTemp(gX.getRows(), gX.getCols());
  CMatrix diagKslash(X.getRows(), 1);
  for(size_t i=0; i<components.size(); i++)
  {
    CTensorKern* slashKern = new CTensorKern(*this, i);
    components[i]->getDiagGradX(gXTemp, X);
    slashKern->diagCompute(diagKslash, X);
    for(unsigned int k=0; k<getInputDim(); k++)
    {
      gXTemp.dotMultiplyColCol(k, diagKslash, 0);
    }
    gX.add(gXTemp);
  }
}
void CTensorKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==getNumParams());
  unsigned int start = 0;
  unsigned int end = 0;
  CMatrix tempCovGrad(covGrad.getRows(), covGrad.getCols());
  for(size_t i=0; i<components.size(); i++)
  {
    CTensorKern* slashKern = new CTensorKern(*this, i);
    slashKern->compute(tempCovGrad, X, X2);
    tempCovGrad.dotMultiply(covGrad);
    end = start+components[i]->getNumParams()-1;
    CMatrix subg(1, components[i]->getNumParams());
    components[i]->getGradParams(subg, X, X2, tempCovGrad, regularise);
    g.setMatrix(0, start, subg);
    start = end+1;
  }
}
void CTensorKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& covGrad, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==getNumParams());
  unsigned int start = 0;
  unsigned int end = 0;
  CMatrix tempCovGrad(covGrad.getRows(), covGrad.getCols());
  for(size_t i=0; i<components.size(); i++)
  {
    CTensorKern* slashKern = new CTensorKern(*this, i);
    slashKern->compute(tempCovGrad, X);
    tempCovGrad.dotMultiply(covGrad);
    end = start+components[i]->getNumParams()-1;
    CMatrix subg(1, components[i]->getNumParams());
    tempCovGrad.setSymmetric(true);
    components[i]->getGradParams(subg, X, tempCovGrad, regularise);
    g.setMatrix(0, start, subg);
    start = end+1;
  }
}
double CTensorKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad) const
{
  
  BOUNDCHECK(index<nParams);
  unsigned int start=0;
  unsigned int end=0;
  CMatrix tempCovGrad(covGrad.getRows(), covGrad.getCols());
  for(size_t i=0; i<components.size(); i++)
  {
    end = start+components[i]->getNumParams()-1;
    if(index<end)
    {
      CTensorKern* slashKern = new CTensorKern(*this, i);
      slashKern->compute(tempCovGrad, X, X2);
      tempCovGrad.dotMultiply(covGrad);
      return components[i]->getGradParam(index-start, X, X2, covGrad);
      start = end + 1;
    }
  }
  return -1;
}
double CTensorKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& covGrad) const
{
  
  BOUNDCHECK(index<nParams);
  unsigned int start=0;
  unsigned int end=0;
  CMatrix tempCovGrad(covGrad.getRows(), covGrad.getCols());
  for(unsigned int i=0; i<components.size(); i++)
  {
    end = start+components[i]->getNumParams()-1;
    if(index<end)
    {
      CTensorKern* slashKern = new CTensorKern(*this, i);
      slashKern->compute(tempCovGrad, X);
      tempCovGrad.dotMultiply(covGrad);
      tempCovGrad.setSymmetric(true);

      return components[i]->getGradParam(index-start, X, covGrad);
      start = end + 1;
    }
  }
  return -1;
}
// the white noise kernel.
CWhiteKern::CWhiteKern() : CKern()
{
  _init();
}
CWhiteKern::CWhiteKern(unsigned int inDim) : CKern(inDim)
{
  _init();
  setInputDim(inDim);
}
CWhiteKern::CWhiteKern(const CMatrix& X) : CKern(X)
{
  _init();
  setInputDim(X.getCols());
}  
CWhiteKern::CWhiteKern(const CWhiteKern& kern) : CKern(kern)
{
  _init();
  setInputDim(kern.getInputDim());
  variance = kern.variance;
}
// Class destructor
CWhiteKern::~CWhiteKern()
{
}
double CWhiteKern::getVariance() const
{
  return variance;
}
void CWhiteKern::_init()
{
  nParams = 1;
  setType("white");
  setName("white noise");
  setParamName("variance", 0);
  addTransform(CTransform::defaultPositive(), 0);
  setStationary(true);
}
void CWhiteKern::setInitParam()
{
  variance = exp(-2.0);  
}

inline double CWhiteKern::diagComputeElement(const CMatrix& X, unsigned int index) const
{
  return variance;
}
void CWhiteKern::diagCompute(CMatrix& d, const CMatrix& X) const
{
  DIMENSIONMATCH(d.getCols()==1);
  DIMENSIONMATCH(X.rowsMatch(d));
  d.setVals(variance);
}
void CWhiteKern::setParam(double val, unsigned int paramNo)
{
  DIMENSIONMATCH(paramNo==0);
  switch(paramNo)
  {
  case 0:
    variance = val;
    break;
  default:
    throw ndlexceptions::Error("Requested parameter doesn't exist.");
  }
}
// Parameters are kernel parameters
double CWhiteKern::getParam(unsigned int paramNo) const
{
  BOUNDCHECK(paramNo==0);
  switch(paramNo)
  {
  case 0:
    return variance;
    break;
  default:
    throw ndlexceptions::Error("Requested parameter doesn't exist.");
  }
}
void CWhiteKern::getGradX(CMatrix& gX, const CMatrix& X, unsigned int row, const CMatrix& X2,  bool addG) const
{
  DIMENSIONMATCH(gX.getRows() == X2.getRows());
  BOUNDCHECK(row < X.getRows());
  DIMENSIONMATCH(X.getCols()==X2.getCols());
  if(!addG)
  {
    gX.zeros();
  }
}
void CWhiteKern::getDiagGradX(CMatrix& gX, const CMatrix& X, bool addG) const
{
  DIMENSIONMATCH(gX.dimensionsMatch(X));
  if(!addG)
    gX.zeros();
}
double CWhiteKern::getWhite() const
{
  return variance;
}

inline double CWhiteKern::computeElement(const CMatrix& X1, unsigned int index1,
					 const CMatrix& X2, unsigned int index2) const
{
  return 0.0;
}

void CWhiteKern::compute(CMatrix& K, const CMatrix& X) const
{
  DIMENSIONMATCH(K.rowsMatch(X));
  MATRIXPROPERTIES(K.isSquare());
  K.zeros();
  for(unsigned int i=0; i<K.getRows(); i++)
    K.setVal(variance, i, i);
  K.setSymmetric(true);
}

void CWhiteKern::compute(CMatrix& K, const CMatrix& X, const CMatrix& X2) const
{
  DIMENSIONMATCH(K.rowsMatch(X));
  DIMENSIONMATCH(K.getCols()==X2.getRows());
  K.zeros();
}
void CWhiteKern::compute(CMatrix& K, const CMatrix& X, const CMatrix& X2, unsigned int row) const
{
  DIMENSIONMATCH(K.rowsMatch(X));
  DIMENSIONMATCH(K.getCols()==1);
  K.zeros();
}
double CWhiteKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad) const
{
  BOUNDCHECK(index==0);
  return 0.0;
}
double CWhiteKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& covGrad) const
{
  BOUNDCHECK(index==0);
  return trace(covGrad);
}

// the white noise kernel.
CWhitefixedKern::CWhitefixedKern() : CKern()
{
  _init();
}
CWhitefixedKern::CWhitefixedKern(unsigned int inDim) : CKern(inDim)
{
  _init();
  setInputDim(inDim);
}
CWhitefixedKern::CWhitefixedKern(const CMatrix& X) : CKern(X)
{
  _init();
  setInputDim(X.getCols());
}  
CWhitefixedKern::CWhitefixedKern(const CWhitefixedKern& kern) : CKern(kern)
{
  _init();
  setInputDim(kern.getInputDim());
  variance = kern.variance;
}
// Class destructor
CWhitefixedKern::~CWhitefixedKern()
{
}
ostream& CWhitefixedKern::display(ostream& os) const
{
  os << getName() << " kernel:" << endl;
  os << "variance: " << variance << endl;
  return os;
}

void CWhitefixedKern::writeParamsToStream(ostream& out) const
{
  writeToStream(out, "baseType", getBaseType());
  writeToStream(out, "type", getType());
  writeToStream(out, "inputDim", getInputDim());
  writeToStream(out, "numParams", getNumParams());
  writeToStream(out, "variance", variance);
}
void CWhitefixedKern::readParamsFromStream(istream& in)
{
//   string tbaseType = getBaseTypeStream(in);
//   if(tbaseType != getBaseType())
//     throw ndlexceptions::StreamFormatError("baseType", "Error mismatch between saved base type, " + tbaseType + ", and Class base type, " + getType() + ".");
//   string ttype = getTypeStream(in);
//   if(ttype != getType())
//     throw ndlexceptions::StreamFormatError("type", "Error mismatch between saved type, " + ttype + ", and Class type, " + getType() + ".");
  setInputDim(readIntFromStream(in, "inputDim"));
  unsigned int numParams=readIntFromStream(in, "numParams");
  setVariance(readDoubleFromStream(in, "variance"));
}

double CWhitefixedKern::getVariance() const
{
  return variance;
}
void CWhitefixedKern::_init()
{
  nParams = 0;
  setType("whitefixed");
  setName("fixed white noise");
  setStationary(true);
}
void CWhitefixedKern::setInitParam()
{
  variance = exp(-2.0);  
}

inline double CWhitefixedKern::diagComputeElement(const CMatrix& X, unsigned int index) const
{
  return variance;
}
void CWhitefixedKern::diagCompute(CMatrix& d, const CMatrix& X) const
{
  DIMENSIONMATCH(d.getCols()==1);
  DIMENSIONMATCH(X.rowsMatch(d));
  d.setVals(variance);
}
void CWhitefixedKern::setParam(double val, unsigned int paramNo)
{
  throw ndlexceptions::Error("Requested parameter doesn't exist.");
}
// Parameters are kernel parameters
double CWhitefixedKern::getParam(unsigned int paramNo) const
{
  throw ndlexceptions::Error("Requested parameter doesn't exist.");
}
void CWhitefixedKern::getGradX(CMatrix& gX, const CMatrix& X, unsigned int row, const CMatrix& X2,  bool addG) const
{
  DIMENSIONMATCH(gX.getRows() == X2.getRows());
  BOUNDCHECK(row < X.getRows());
  DIMENSIONMATCH(X.getCols()==X2.getCols());
  if(!addG)
  {
    gX.zeros();
  }
}
void CWhitefixedKern::getDiagGradX(CMatrix& gX, const CMatrix& X, bool addG) const
{
  DIMENSIONMATCH(gX.dimensionsMatch(X));
  if(!addG)
    gX.zeros();
}
double CWhitefixedKern::getWhite() const
{
  return variance;
}

inline double CWhitefixedKern::computeElement(const CMatrix& X1, unsigned int index1,
					 const CMatrix& X2, unsigned int index2) const
{
  return 0.0;
}

void CWhitefixedKern::compute(CMatrix& K, const CMatrix& X) const
{
  DIMENSIONMATCH(K.rowsMatch(X));
  MATRIXPROPERTIES(K.isSquare());
  K.zeros();
  for(unsigned int i=0; i<K.getRows(); i++)
    K.setVal(variance, i, i);
  K.setSymmetric(true);
}

void CWhitefixedKern::compute(CMatrix& K, const CMatrix& X, const CMatrix& X2) const
{
  DIMENSIONMATCH(K.rowsMatch(X));
  DIMENSIONMATCH(K.getCols()==X2.getRows());
  K.zeros();
}
void CWhitefixedKern::compute(CMatrix& K, const CMatrix& X, const CMatrix& X2, unsigned int row) const
{
  DIMENSIONMATCH(K.rowsMatch(X));
  DIMENSIONMATCH(K.getCols()==1);
  K.zeros();
}
double CWhitefixedKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad) const
{
  throw ndlexceptions::Error("Requested parameter doesn't exist.");

}
double CWhitefixedKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& covGrad) const
{
  throw ndlexceptions::Error("Requested parameter doesn't exist.");
}

// the bias kernel.
CBiasKern::CBiasKern() : CKern()
{
  _init();
}
CBiasKern::CBiasKern(unsigned int inDim) : CKern(inDim)
{
  _init();
  setInputDim(inDim);
}
CBiasKern::CBiasKern(const CMatrix& X) : CKern(X)
{
  _init();
  setInputDim(X.getCols());
}  
CBiasKern::CBiasKern(const CBiasKern& kern) : CKern(kern)
{
  _init();
  setInputDim(kern.getInputDim());
  variance = kern.variance;
  
}
// Class destructor
CBiasKern::~CBiasKern()
{
}
double CBiasKern::getVariance() const
{
  return variance;
}
void CBiasKern::_init()
{
  nParams = 1;
  setType("bias");
  setName("bias");
  setParamName("variance", 0);
  addTransform(CTransform::defaultPositive(), 0);
  setStationary(true);
}

void CBiasKern::setInitParam()
{
  variance = exp(-2.0);  
}

double CBiasKern::diagComputeElement(const CMatrix& X, unsigned int index) const
{
  return variance;
}
void CBiasKern::diagCompute(CMatrix& d, const CMatrix& X) const
{
  DIMENSIONMATCH(d.getCols()==1);
  DIMENSIONMATCH(X.rowsMatch(d));
  d.setVals(variance);
}
// Parameters are kernel parameters
void CBiasKern::setParam(double val, unsigned int paramNo)
{
  BOUNDCHECK(paramNo==0);
  switch(paramNo)
  {
  case 0:
    variance = val;
    break;
  default:
    throw ndlexceptions::Error("Requested parameter doesn't exist.");
  }
}
double CBiasKern::getParam(unsigned int paramNo) const
{
  BOUNDCHECK(paramNo==0);
  switch(paramNo)
  {
  case 0:
    return variance;
    break;
  default:
    throw ndlexceptions::Error("Requested parameter doesn't exist.");
  }
}
void CBiasKern::getGradX(CMatrix& gX, const CMatrix& X, unsigned int row, const CMatrix& X2, bool addG) const
{
  DIMENSIONMATCH(gX.getRows() == X2.getRows());
  BOUNDCHECK(row < X.getRows());
  DIMENSIONMATCH(X.getCols()==X2.getCols());
  if(!addG)
  {
    gX.zeros();
  }
}
void CBiasKern::getDiagGradX(CMatrix& gX, const CMatrix& X, bool addG) const
{
  DIMENSIONMATCH(gX.dimensionsMatch(X));
  if(!addG)
    gX.zeros();
}
double CBiasKern::getWhite() const
{
  return 0.0;
}

inline double CBiasKern::computeElement(const CMatrix& X1, unsigned int index1, 
					const CMatrix& X2, unsigned int index2) const
{
  return variance;
}

void CBiasKern::compute(CMatrix& K, const CMatrix& X) const
{
  DIMENSIONMATCH(K.rowsMatch(X));
  MATRIXPROPERTIES(K.isSquare());
  K.setVals(variance);
  K.setSymmetric(true);
}

void CBiasKern::compute(CMatrix& K, const CMatrix& X, const CMatrix& X2) const
{
  DIMENSIONMATCH(K.rowsMatch(X));
  DIMENSIONMATCH(K.getCols()==X2.getRows());
  K.setVals(variance);
}
void CBiasKern::compute(CMatrix& K, const CMatrix& X, const CMatrix& X2, unsigned int row) const
{
  DIMENSIONMATCH(K.rowsMatch(X));
  DIMENSIONMATCH(K.getCols()==1);
  K.setVals(variance);
}
double CBiasKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad) const 
{
  BOUNDCHECK(index==0);
  return covGrad.sum();
}
double CBiasKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& covGrad) const 
{
  BOUNDCHECK(index==0);
  return covGrad.sum();
}
// the RBF kernel.
CRbfKern::CRbfKern() : CKern()
{
  _init();
}
CRbfKern::CRbfKern(unsigned int inDim) : CKern(inDim)
{
  _init();
  setInputDim(inDim);
}
CRbfKern::CRbfKern(const CMatrix& X) : CKern(X)
{
  _init();
  setInputDim(X.getCols());
}  
CRbfKern::CRbfKern(const CRbfKern& kern) : CKern(kern)
{
  _init();
  setInputDim(kern.getInputDim());
  variance = kern.variance;
  inverseWidth = kern.inverseWidth;
}
// Class destructor
CRbfKern::~CRbfKern()
{
}
double CRbfKern::getVariance() const
{
  return variance;
}
void CRbfKern::_init()
{
  nParams = 2;
  setType("rbf");
  setName("RBF");
  setParamName("inverseWidth", 0);
  addTransform(CTransform::defaultPositive(), 0);
  setParamName("variance", 1);
  addTransform(CTransform::defaultPositive(), 1);
  setStationary(true);
}
void CRbfKern::setInitParam()
{
  inverseWidth = 1.0;
  variance = 1.0;
}

inline double CRbfKern::diagComputeElement(const CMatrix& X, unsigned int index) const
{
  return variance;
}
void CRbfKern::diagCompute(CMatrix& d, const CMatrix& X) const
{
  DIMENSIONMATCH(d.getCols()==1);
  DIMENSIONMATCH(X.rowsMatch(d));
  d.setVals(variance);
}
// Parameters are kernel parameters
void CRbfKern::setParam(double val, unsigned int paramNo)
{
  BOUNDCHECK(paramNo < nParams);
  switch(paramNo)
  {
  case 0:
    inverseWidth = val;
    break;
  case 1:
    variance = val;
    break;
  default:
    throw ndlexceptions::Error("Requested parameter doesn't exist.");
  }
}
double CRbfKern::getParam(unsigned int paramNo) const
{
  BOUNDCHECK(paramNo < nParams);
  switch(paramNo)
  {
  case 0:
    return inverseWidth;
    break;
  case 1:
    return variance;
    break;
  default:
    throw ndlexceptions::Error("Requested parameter doesn't exist.");
  }
}
void CRbfKern::getGradX(CMatrix& gX, const CMatrix& X, unsigned int row, const CMatrix& X2, bool addG) const
{
  DIMENSIONMATCH(gX.getRows() == X2.getRows());
  BOUNDCHECK(row < X.getRows());
  DIMENSIONMATCH(X.getCols()==X2.getCols());
  double wi2 = 0.5*inverseWidth;
  double pf = variance*inverseWidth;
  DIMENSIONMATCH(gX.getCols()==X2.getCols());
  for(unsigned int k=0; k<X2.getRows(); k++)
  {
    double n2 = X.dist2Row(row, X2, k);
    for(unsigned int j=0; j<X2.getCols(); j++)
    {
      double val = pf*(X2.getVal(k, j)-X.getVal(row, j))*exp(-n2*wi2);
      if(addG)
	gX.addVal(val, k, j);
      else 
	gX.setVal(val, k, j);
    }
  }
}
void CRbfKern::getDiagGradX(CMatrix& gX, const CMatrix& X, bool addG) const
{
  DIMENSIONMATCH(gX.dimensionsMatch(X));
  if(!addG)
    gX.zeros();
}
double CRbfKern::getWhite() const
{
  return 0.0;
}

double CRbfKern::computeElement(const CMatrix& X1, unsigned int index1, 
				const CMatrix& X2, unsigned int index2) const
{
  double k = X1.dist2Row(index1, X2, index2);
  k = 0.5*k*inverseWidth;
  k = variance*exp(-k);
  return k;
}

void CRbfKern::updateX(const CMatrix& X)
{
  setUpdateXused(true);
  Xdists.resize(X.getRows(),X.getRows());
  double halfInverseWidth=0.5*inverseWidth;
  unsigned int nrows = X.getRows();
  for(unsigned int j=0; j<nrows; j++)
  {
    Xdists.setVal(0,j,j);
    for(unsigned int i=0; i<j; i++)
    {
      double dist2 = X.dist2Row(i, X, j);
      Xdists.setVal(dist2,i,j);
      Xdists.setVal(exp(-dist2*halfInverseWidth),j,i);
    }
  }
}
void CRbfKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==nParams);
  DIMENSIONMATCH(X.getRows()==covGrad.getRows());
  DIMENSIONMATCH(X2.getRows()==covGrad.getCols());
  double g1=0.0;
  double g2=0.0;
  double halfInverseWidth=0.5*inverseWidth;
  double halfVariance=0.5*variance;

  unsigned int nrows = X.getRows();
  for(unsigned int j=0; j<nrows; j++)
  {
    for(unsigned int i=0; i<X2.getRows(); i++)
    {
      double k = 0;
      double dist2 = 0;
      dist2 = X2.dist2Row(i, X, j);
      k = exp(-dist2*halfInverseWidth);
      double kcg_ij = k*covGrad.getVal(j,i);
      g1 -= halfVariance*dist2*kcg_ij; // dk()/dgamma in SBIK paper
      g2 += kcg_ij;                    // dk()/dalpha in SBIK paper
    }
  }
  g.setVal(g1, 0);
  g.setVal(g2, 1);
  if(regularise)
    addPriorGrad(g);
}

void CRbfKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& covGrad, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==nParams);
  MATRIXPROPERTIES(covGrad.isSymmetric());
  double g1=0.0;
  double g2=0.0;
  double halfInverseWidth=0.5*inverseWidth;
  double halfVariance=0.5*variance;

  unsigned int nrows = X.getRows();
  for(unsigned int j=0; j<nrows; j++)
  {
    g2 += covGrad.getVal(j,j);
    for(unsigned int i=0; i<j; i++)
    {
      double k = 0;
      double dist2 = 0;
      if(isUpdateXused()) // WVB's mod for precomputing parts of the kernel.
      {
	dist2 = Xdists.getVal(i,j);
	k = Xdists.getVal(j,i);
      }
      else
      {
	dist2 = X.dist2Row(i, X, j);
	k = exp(-dist2*halfInverseWidth);
      }
      double kcg_ij = k*covGrad.getVal(i,j);
      g1 -= 2.0*halfVariance*dist2*kcg_ij; // dk()/dgamma in SBIK paper
      g2 += 2.0*kcg_ij;                    // dk()/dalpha in SBIK paper
    }
  }
  g.setVal(g1, 0);
  g.setVal(g2, 1);
  if(regularise)
    addPriorGrad(g);
}

double CRbfKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad) const
{
  
  BOUNDCHECK(index<nParams);
  throw ndlexceptions::NotImplementedError( "Error getGradParam is not currently implemented for CRbfKern");
  
}
double CRbfKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& covGrad) const
{
  
  BOUNDCHECK(index<nParams);
  throw ndlexceptions::NotImplementedError( "Error getGradParam is not currently implemented for CRbfKern");
  
}

// the Rational Quadratic kernel.
CRatQuadKern::CRatQuadKern() : CKern()
{
  _init();
}
CRatQuadKern::CRatQuadKern(unsigned int inDim) : CKern(inDim)
{
  _init();
  setInputDim(inDim);
}
CRatQuadKern::CRatQuadKern(const CMatrix& X) : CKern(X)
{
  _init();
  setInputDim(X.getCols());
}  
CRatQuadKern::CRatQuadKern(const CRatQuadKern& kern) : CKern(kern)
{
  _init();
  setInputDim(kern.getInputDim());
  variance = kern.variance;
  alpha = kern.alpha;
  lengthScale = kern.lengthScale;
}
// Class destructor
CRatQuadKern::~CRatQuadKern()
{
}
double CRatQuadKern::getVariance() const
{
  return variance;
}
void CRatQuadKern::_init()
{
  nParams = 3;
  setType("ratquad");
  setName("Rational Quadratic");
  setParamName("alpha", 0);
  addTransform(CTransform::defaultPositive(), 0);
  setParamName("lengthScale", 1);
  addTransform(CTransform::defaultPositive(), 1);
  setParamName("variance", 2);
  addTransform(CTransform::defaultPositive(), 2);
  setStationary(true);
}
void CRatQuadKern::setInitParam()
{
  alpha = 1.0;
  lengthScale = 1.0;
  variance = 1.0;
}

inline double CRatQuadKern::diagComputeElement(const CMatrix& X, unsigned int index) const
{
  return variance;
}
void CRatQuadKern::diagCompute(CMatrix& d, const CMatrix& X) const
{
  DIMENSIONMATCH(d.getCols()==1);
  DIMENSIONMATCH(X.rowsMatch(d));
  d.setVals(variance);
}
// Parameters are kernel parameters
void CRatQuadKern::setParam(double val, unsigned int paramNo)
{
  BOUNDCHECK(paramNo < nParams);
  switch(paramNo)
  {
  case 0:
    alpha = val;
    break;
  case 1:
    lengthScale = val;
    break;
  case 2:
    variance = val;
    break;
  default:
    throw ndlexceptions::Error("Requested parameter doesn't exist.");
  }
}
double CRatQuadKern::getParam(unsigned int paramNo) const
{
  BOUNDCHECK(paramNo < nParams);
  switch(paramNo)
  {
  case 0:
    return alpha;
    break;
  case 1:
    return lengthScale;
    break;
  case 2:
    return variance;
    break;
  default:
    throw ndlexceptions::Error("Requested parameter doesn't exist.");
  }
}

void CRatQuadKern::getGradX(CMatrix& gX, const CMatrix& X, unsigned int row, const CMatrix& X2, bool addG) const
{
  DIMENSIONMATCH(gX.getRows() == X2.getRows());
  BOUNDCHECK(row < X.getRows());
  DIMENSIONMATCH(X.getCols()==X2.getCols());
  double wi2 = 0.5/(lengthScale*lengthScale*alpha);
  for(unsigned int k=0; k<X2.getRows(); k++)
  {
    double n2 = X.dist2Row(row, X2, k);
    double ratquadPart = variance*pow((1+n2*wi2),-(alpha+1))/(lengthScale*lengthScale);
    for(unsigned int j=0; j<X2.getCols(); j++)
    {
      double val = ratquadPart*(X2.getVal(k, j)-X.getVal(row, j));
      if(addG)
	gX.addVal(val, k, j);
      else 
	gX.setVal(val, k, j);
    }
  }
}
void CRatQuadKern::getDiagGradX(CMatrix& gX, const CMatrix& X, bool addG) const
{
  DIMENSIONMATCH(gX.dimensionsMatch(X));
  if(!addG)
    gX.zeros();
}
double CRatQuadKern::getWhite() const
{
  return 0.0;
}

double CRatQuadKern::computeElement(const CMatrix& X1, unsigned int index1, 
				    const CMatrix& X2, unsigned int index2) const
{
  double k = X1.dist2Row(index1, X2, index2);
  k *= 0.5/(lengthScale*lengthScale*alpha);
  k = variance*pow((1+k), -alpha);
  return k;
}

void CRatQuadKern::updateX(const CMatrix& X)
{
  setUpdateXused(true);
  Xdists.resize(X.getRows(),X.getRows());
  double wi2=0.5/(lengthScale*lengthScale*alpha);
  unsigned int nrows = X.getRows();
  for(unsigned int j=0; j<nrows; j++)
  {
    Xdists.setVal(0,j,j);
    for(unsigned int i=0; i<j; i++)
    {
      double dist2 = X.dist2Row(i, X, j);
      Xdists.setVal(dist2,i,j);
      Xdists.setVal(pow((1+dist2*wi2), -alpha),j,i);
    }
  }
}
void CRatQuadKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==nParams);
  DIMENSIONMATCH(X.getRows()==covGrad.getRows());
  DIMENSIONMATCH(X2.getRows()==covGrad.getCols());
  double g1=0.0;
  double g2=0.0;
  double g3=0.0;
  double wi2=0.5/(lengthScale*lengthScale*alpha);

  unsigned int nrows = X.getRows();
  for(unsigned int j=0; j<nrows; j++)
  {
    for(unsigned int i=0; i<X2.getRows(); i++)
    {
      double dist2 = X2.dist2Row(i, X, j);
      double baseVal = (1+dist2*wi2);
      double kbase = pow(baseVal, -alpha);
      double kbase2 = -alpha*kbase/baseVal;
      double cgVal = covGrad.getVal(j,i);
      double kcg_ij = kbase*cgVal;
      g1 -= cgVal*(dist2*wi2/alpha*kbase2 + log(baseVal)*kbase); 
      g2 -= cgVal*(dist2*kbase2);
      g3 += kcg_ij;                   
    }
  }
  g.setVal(g1*variance, 0);
  g.setVal(g2*variance*2*wi2/lengthScale, 1);
  g.setVal(g3, 2);
  if(regularise)
    addPriorGrad(g);
}

void CRatQuadKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& covGrad, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==nParams);
  MATRIXPROPERTIES(covGrad.isSymmetric());
  double g1=0.0;
  double g2=0.0;
  double g3=0.0;
  double wi2=0.5/(lengthScale*lengthScale*alpha);
  
  unsigned int nrows = X.getRows();
  for(unsigned int j=0; j<nrows; j++)
  {
    double baseVal = 0.0;
    double kbase = 0.0;
    double kbase2 = 0.0;    
    double dist2 = 0.0;
    for(unsigned int i=0; i<j; i++)
    {
      if(isUpdateXused()) // WVB's mod for precomputing parts of the kernel.
      {
	dist2 = Xdists.getVal(i,j);
	baseVal = (1+dist2*wi2);
	kbase = Xdists.getVal(j,i);
      }
      else
      {
	dist2 = X.dist2Row(i, X, j);
	baseVal = (1+dist2*wi2);
	kbase = pow(baseVal, -alpha);
      }
      kbase2 = -alpha*kbase/baseVal;
      double cgVal = covGrad.getVal(j,i);
      double kcg_ij = kbase*cgVal;
      g1 -= cgVal*(dist2*wi2/alpha*kbase2 + log(baseVal)*kbase); 
      g2 -= cgVal*(dist2*kbase2);
      g3 += kcg_ij;                   
    }
    g3+=0.5*covGrad.getVal(j, j);
  }
  g.setVal(2.0*g1*variance, 0);
  g.setVal(4.0*g2*variance*wi2/lengthScale, 1);
  g.setVal(2.0*g3, 2);
  if(regularise)
    addPriorGrad(g);
}

double CRatQuadKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad) const
{
  
  BOUNDCHECK(index<nParams);
  throw ndlexceptions::NotImplementedError("getGradParam is not currently implemented for CRatQuadKern.");
}
double CRatQuadKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& covGrad) const
{
  
  BOUNDCHECK(index<nParams);
  throw ndlexceptions::NotImplementedError("getGradParam is not currently implemented for CRatQuadKern.");
}


// the Matern 3/2 kernel.
CMatern32Kern::CMatern32Kern() : CKern()
{
  _init();
}
CMatern32Kern::CMatern32Kern(unsigned int inDim) : CKern(inDim)
{
  _init();
  setInputDim(inDim);
}
CMatern32Kern::CMatern32Kern(const CMatrix& X) : CKern(X)
{
  _init();
  setInputDim(X.getCols());
}  
CMatern32Kern::CMatern32Kern(const CMatern32Kern& kern) : CKern(kern)
{
  _init();
  setInputDim(kern.getInputDim());
  variance = kern.variance;
  lengthScale = kern.lengthScale;
}
// Class destructor
CMatern32Kern::~CMatern32Kern()
{
}
double CMatern32Kern::getVariance() const
{
  return variance;
}
void CMatern32Kern::_init()
{
  nParams = 2;
  setType("matern32");
  setName("Matern dof=3/2");
  setParamName("lengthScale", 0);
  addTransform(CTransform::defaultPositive(), 0);
  setParamName("variance", 1);
  addTransform(CTransform::defaultPositive(), 1);
  setStationary(true);
}
void CMatern32Kern::setInitParam()
{
  lengthScale = 1.0;
  variance = 1.0;
}

inline double CMatern32Kern::diagComputeElement(const CMatrix& X, unsigned int index) const
{
  return variance;
}
void CMatern32Kern::diagCompute(CMatrix& d, const CMatrix& X) const
{
  DIMENSIONMATCH(d.getCols()==1);
  DIMENSIONMATCH(X.rowsMatch(d));
  d.setVals(variance);
}
// Parameters are kernel parameters
void CMatern32Kern::setParam(double val, unsigned int paramNo)
{
  BOUNDCHECK(paramNo < nParams);
  switch(paramNo)
  {
  case 0:
    lengthScale = val;
    break;
  case 1:
    variance = val;
    break;
  default:
    throw ndlexceptions::Error("Requested parameter doesn't exist.");
  }
}
double CMatern32Kern::getParam(unsigned int paramNo) const
{
  BOUNDCHECK(paramNo < nParams);
  switch(paramNo)
  {
  case 0:
    return lengthScale;
    break;
  case 1:
    return variance;
    break;
  default:
    throw ndlexceptions::Error("Requested parameter doesn't exist.");
  }
}

void CMatern32Kern::getGradX(CMatrix& gX, const CMatrix& X, unsigned int row, const CMatrix& X2, bool addG) const
{
  DIMENSIONMATCH(gX.getRows() == X2.getRows());
  BOUNDCHECK(row < X.getRows());
  DIMENSIONMATCH(X.getCols()==X2.getCols());
  double wi2= 3.0/(lengthScale*lengthScale);
  for(unsigned int k=0; k<X2.getRows(); k++)
  {
    double n2 = X.dist2Row(row, X2, k);
    double sqrtn2wi2 = sqrt(n2*wi2);
    double expmSqrtn2wi2 = exp(-sqrtn2wi2);
    double K = variance*(1+sqrtn2wi2)*expmSqrtn2wi2;
    for(unsigned int j=0; j<X2.getCols(); j++)
    {
      double ratio;
      if(sqrtn2wi2 != 0)
	ratio = (X2.getVal(k, j)-X.getVal(row, j))/sqrtn2wi2;
      else
	ratio = 1.0;
      double val = wi2*ratio*(K-variance*expmSqrtn2wi2);
      if(addG)
	gX.addVal(val, k, j);
      else 
	gX.setVal(val, k, j);
    }
  }
}
void CMatern32Kern::getDiagGradX(CMatrix& gX, const CMatrix& X, bool addG) const
{
  DIMENSIONMATCH(gX.dimensionsMatch(X));
  if(!addG)
    gX.zeros();
}
double CMatern32Kern::getWhite() const
{
  return 0.0;
}

double CMatern32Kern::computeElement(const CMatrix& X1, unsigned int index1, 
				     const CMatrix& X2, unsigned int index2) const
{
  double k = X1.dist2Row(index1, X2, index2);
  double wi2 = (3.0/(lengthScale*lengthScale));
  double sqrtn2wi2 = sqrt(k*wi2);
  k = variance*(1+sqrtn2wi2)*exp(-sqrtn2wi2);
  return k;
}

void CMatern32Kern::updateX(const CMatrix& X)
{
  setUpdateXused(true);
  Xdists.resize(X.getRows(),X.getRows());
  double wi2=3.0/(lengthScale*lengthScale);
  unsigned int nrows = X.getRows();
  for(unsigned int j=0; j<nrows; j++)
  {
    Xdists.setVal(0,j,j);
    for(unsigned int i=0; i<j; i++)
    {
      double n2 = X.dist2Row(i, X, j);
      Xdists.setVal(n2,i,j);
      Xdists.setVal(sqrt(n2*wi2),j,i);
    }
  }
}
void CMatern32Kern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==nParams);
  DIMENSIONMATCH(X.getRows()==covGrad.getRows());
  DIMENSIONMATCH(X2.getRows()==covGrad.getCols());
  double g1=0.0;
  double g2=0.0;
  double wi2=3.0/(lengthScale*lengthScale);
  unsigned int nrows = X.getRows();
  for(unsigned int j=0; j<nrows; j++)
  {
    for(unsigned int i=0; i<X2.getRows(); i++)
    {
      double cgVal = covGrad.getVal(j,i);
      double n2 = X2.dist2Row(i, X, j);
      double sqrtn2wi2 = sqrt(n2*wi2);
      double expmsqrtn2wi2 = exp(-sqrtn2wi2);
      double k = (1+sqrtn2wi2)*expmsqrtn2wi2;
      double ratio;
      if(sqrtn2wi2==0.0)
	ratio = 1.0;
      else
	ratio = n2/sqrtn2wi2;
      g1 += cgVal*wi2/lengthScale*ratio*(k-expmsqrtn2wi2);
      g2 += cgVal*k;
    }
  }
  g.setVal(g1*variance, 0);
  g.setVal(g2, 1);
  if(regularise)
    addPriorGrad(g);
}

void CMatern32Kern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& covGrad, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==nParams);
  MATRIXPROPERTIES(covGrad.isSymmetric());
  double g1=0.0;
  double g2=0.0;
  double wi2=3.0/(lengthScale*lengthScale);
  double n2 = 0.0;
  double sqrtwi2n2 = 0.0;
  double ratio;

  unsigned int nrows = X.getRows();
  for(unsigned int j=0; j<nrows; j++)
  {
    for(unsigned int i=0; i<j; i++)
    {
      if(isUpdateXused()) // WVB's mod for precomputing parts of the kernel.
      {
	n2 = Xdists.getVal(i,j);
	sqrtwi2n2 = Xdists.getVal(j,i);
      }
      else
      {
	n2 = X.dist2Row(i, X, j);
	sqrtwi2n2 = sqrt(n2*wi2);
      }
      if(sqrtwi2n2==0.0)
	ratio = 1.0;
      else
	ratio = n2/sqrtwi2n2;
      double expmsqrtwi2n2 = exp(-sqrt(n2*wi2));
      double cgVal = covGrad.getVal(j,i);
      double k = (1+sqrtwi2n2)*expmsqrtwi2n2;
      g1 += cgVal*wi2/lengthScale*ratio*(k-expmsqrtwi2n2);
      g2 += cgVal*k;
    }
    g2 += 0.5*covGrad.getVal(j, j);
  }
  g.setVal(2.0*g1*variance, 0);
  g.setVal(2.0*g2, 1);
  if(regularise)
    addPriorGrad(g);
}

double CMatern32Kern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad) const
{
  
  BOUNDCHECK(index<nParams);
  throw ndlexceptions::NotImplementedError("Error getGradParam is not currently implemented for CMatern32Kern");
  
}
double CMatern32Kern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& covGrad) const
{
  
  BOUNDCHECK(index<nParams);
  throw ndlexceptions::NotImplementedError("Error getGradParam is not currently implemented for CMatern32Kern"); 
}

// the Matern 5/2 kernel.
CMatern52Kern::CMatern52Kern() : CKern()
{
  _init();
}
CMatern52Kern::CMatern52Kern(unsigned int inDim) : CKern(inDim)
{
  _init();
  setInputDim(inDim);
}
CMatern52Kern::CMatern52Kern(const CMatrix& X) : CKern(X)
{
  _init();
  setInputDim(X.getCols());
}  
CMatern52Kern::CMatern52Kern(const CMatern52Kern& kern) : CKern(kern)
{
  _init();
  setInputDim(kern.getInputDim());
  variance = kern.variance;
  lengthScale = kern.lengthScale;
}
// Class destructor
CMatern52Kern::~CMatern52Kern()
{
}
double CMatern52Kern::getVariance() const
{
  return variance;
}
void CMatern52Kern::_init()
{
  nParams = 2;
  setType("matern52");
  setName("Matern dof=5/2");
  setParamName("lengthScale", 0);
  addTransform(CTransform::defaultPositive(), 0);
  setParamName("variance", 1);
  addTransform(CTransform::defaultPositive(), 1);
  setStationary(true);
}
void CMatern52Kern::setInitParam()
{
  lengthScale = 1.0;
  variance = 1.0;
}

inline double CMatern52Kern::diagComputeElement(const CMatrix& X, unsigned int index) const
{
  return variance;
}
void CMatern52Kern::diagCompute(CMatrix& d, const CMatrix& X) const
{
  DIMENSIONMATCH(d.getCols()==1);
  DIMENSIONMATCH(X.rowsMatch(d));
  d.setVals(variance);
}
// Parameters are kernel parameters
void CMatern52Kern::setParam(double val, unsigned int paramNo)
{
  BOUNDCHECK(paramNo < nParams);
  switch(paramNo)
  {
  case 0:
    lengthScale = val;
    break;
  case 1:
    variance = val;
    break;
  default:
    throw ndlexceptions::Error("Requested parameter doesn't exist.");
  }
}
double CMatern52Kern::getParam(unsigned int paramNo) const
{
  BOUNDCHECK(paramNo < nParams);
  switch(paramNo)
  {
  case 0:
    return lengthScale;
    break;
  case 1:
    return variance;
    break;
  default:
    throw ndlexceptions::Error("Requested parameter doesn't exist.");
  }
}
void CMatern52Kern::getGradX(CMatrix& gX, const CMatrix& X, unsigned int row, const CMatrix& X2, bool addG) const
{
  DIMENSIONMATCH(gX.getRows() == X2.getRows());
  BOUNDCHECK(row < X.getRows());
  DIMENSIONMATCH(X.getCols()==X2.getCols());
  double wi2= 5.0/(lengthScale*lengthScale);
  double n2 = 0.0;
  double n2wi2 = 0.0;
  double sqrtn2wi2 = 0.0;
  double K = 0.0;
  double ratio = 0.0;
  double val = 0.0;
  double expmSqrtn2wi2 = 0.0;
  for(unsigned int k=0; k<X2.getRows(); k++)
  {
    n2 = X.dist2Row(row, X2, k);
    n2wi2 = n2*wi2;
    sqrtn2wi2 = sqrt(n2wi2);
    expmSqrtn2wi2 = exp(-sqrtn2wi2);
    K = variance*(1+sqrtn2wi2+n2wi2/3.0)*expmSqrtn2wi2;
    for(unsigned int j=0; j<X2.getCols(); j++)
    {
      if(sqrtn2wi2 != 0.0)
	ratio = (X2.getVal(k, j)-X.getVal(row, j))/sqrtn2wi2;
      else
	ratio = 1.0;
      val = wi2*ratio*(K-variance*expmSqrtn2wi2)-variance*2.0*wi2/3.0*(X2.getVal(k, j)-X.getVal(row, j))*expmSqrtn2wi2;
      if(addG)
	gX.addVal(val, k, j);
      else 
	gX.setVal(val, k, j);
    }
  } 
}
void CMatern52Kern::getDiagGradX(CMatrix& gX, const CMatrix& X, bool addG) const
{
  DIMENSIONMATCH(gX.dimensionsMatch(X));
  if(!addG)
    gX.zeros();
}
double CMatern52Kern::getWhite() const
{
  return 0.0;
}

double CMatern52Kern::computeElement(const CMatrix& X1, unsigned int index1, 
				     const CMatrix& X2, unsigned int index2) const
{
  double k = X1.dist2Row(index1, X2, index2);
  double wi2 = (5.0/(lengthScale*lengthScale));
  double  n2wi2 = k*wi2;
  double sqrtn2wi2 = sqrt(n2wi2);
  k = variance*(1+sqrtn2wi2+n2wi2/3.0)*exp(-sqrtn2wi2);
  return k;
}

void CMatern52Kern::updateX(const CMatrix& X)
{
  setUpdateXused(true);
  Xdists.resize(X.getRows(),X.getRows());
  double wi2=5.0/(lengthScale*lengthScale);
  unsigned int nrows = X.getRows();
  for(unsigned int j=0; j<nrows; j++)
  {
    Xdists.setVal(0,j,j);
    for(unsigned int i=0; i<j; i++)
    {
      double n2 = X.dist2Row(i, X, j);
      Xdists.setVal(n2,i,j);
      Xdists.setVal(sqrt(n2*wi2),j,i);
    }
  }
}
void CMatern52Kern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==nParams);
  DIMENSIONMATCH(X.getRows()==covGrad.getRows());
  DIMENSIONMATCH(X2.getRows()==covGrad.getCols());
  double g1=0.0;
  double g2=0.0;
  double wi2=5.0/(lengthScale*lengthScale);
  double n2wi2 = 0.0;
  double sqrtn2wi2 = 0.0;
  double expmsqrtn2wi2 = 0.0;
  double n2 = 0.0;
  double k = 0.0;
  double cgVal = 0.0;
  double ratio = 0.0;
  unsigned int nrows = X.getRows();
  for(unsigned int j=0; j<nrows; j++)
  {
    for(unsigned int i=0; i<X2.getRows(); i++)
    {
      cgVal = covGrad.getVal(j,i);
      n2 = X2.dist2Row(i, X, j);
      n2wi2 = n2*wi2;
      sqrtn2wi2 = sqrt(n2wi2);
      expmsqrtn2wi2 = exp(-sqrtn2wi2);
      k = (1+sqrtn2wi2+n2wi2/3.0)*expmsqrtn2wi2;
      if(sqrtn2wi2==0.0)
	ratio = 1.0;
      else
	ratio = n2/sqrtn2wi2;
      g1 += cgVal*(wi2/lengthScale*ratio*(k-expmsqrtn2wi2)-2.0*wi2/(3.0*lengthScale)*n2*expmsqrtn2wi2);
      g2 += cgVal*k;
    }
  }
  g.setVal(g1*variance, 0);
  g.setVal(g2, 1);
  if(regularise)
    addPriorGrad(g);
}

void CMatern52Kern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& covGrad, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==nParams);
  MATRIXPROPERTIES(covGrad.isSymmetric());
  double g1=0.0;
  double g2=0.0;
  double wi2=5.0/(lengthScale*lengthScale);
  double n2 = 0.0;
  double n2wi2 = 0.0;
  double sqrtwi2n2 = 0.0;
  double ratio;

  unsigned int nrows = X.getRows();
  for(unsigned int j=0; j<nrows; j++)
  {
    for(unsigned int i=0; i<j; i++)
    {
      if(isUpdateXused()) // WVB's mod for precomputing parts of the kernel.
      {
	n2 = Xdists.getVal(i,j);
	n2wi2 = n2*wi2;
	sqrtwi2n2 = Xdists.getVal(j,i);
      }
      else
      {
	n2 = X.dist2Row(i, X, j);
	n2wi2 = n2*wi2;
	sqrtwi2n2 = sqrt(n2wi2);
      }
      if(sqrtwi2n2==0.0)
	ratio = 1.0;
      else
	ratio = n2/sqrtwi2n2;
      double expmsqrtwi2n2 = exp(-sqrtwi2n2);
      double cgVal = covGrad.getVal(j,i);
      double k = (1+sqrtwi2n2+n2wi2/3.0)*expmsqrtwi2n2;
      g1 += cgVal*(wi2/lengthScale*ratio*(k-expmsqrtwi2n2)-2.0*wi2/(3.0*lengthScale)*n2*expmsqrtwi2n2);
      g2 += cgVal*k;
    }
    g2 += 0.5*covGrad.getVal(j, j);
  }
  g.setVal(2.0*g1*variance, 0);
  g.setVal(2.0*g2, 1);
  if(regularise)
    addPriorGrad(g);
}

double CMatern52Kern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad) const
{
  
  BOUNDCHECK(index<nParams);
  throw ndlexceptions::NotImplementedError("Error getGradParam is not currently implemented for CMatern52Kern");
  
}
double CMatern52Kern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& covGrad) const
{
  
  BOUNDCHECK(index<nParams);
  throw ndlexceptions::NotImplementedError("Error getGradParam is not currently implemented for CMatern52Kern"); 
}


// the Linear kernel.
CLinKern::CLinKern() : CKern()
{
  _init();
}
CLinKern::CLinKern(unsigned int inDim) : CKern(inDim)
{
  _init();
  setInputDim(inDim);
}
CLinKern::CLinKern(const CMatrix& X) : CKern(X)
{
  _init();
  setInputDim(X.getCols());
}  
CLinKern::CLinKern(const CLinKern& kern) : CKern(kern)
{
  _init();
  setInputDim(kern.getInputDim());
  variance = kern.variance;
}
// Class destructor
CLinKern::~CLinKern()
{
}
double CLinKern::getVariance() const
{
  return variance;
}
void CLinKern::_init()
{
  nParams = 1;
  setType("lin");
  setName("linear");
  setParamName("variance", 0);
  addTransform(CTransform::defaultPositive(), 0);
  setStationary(false);
}
void CLinKern::setInitParam()
{
  variance = 1.0;
}

double CLinKern::diagComputeElement(const CMatrix& X, unsigned int index1) const
{
  return variance*X.norm2Row(index1);  
}
// Parameters are kernel parameters
void CLinKern::setParam(double val, unsigned int paramNo)
{
  BOUNDCHECK(paramNo==0);
  switch(paramNo)
  {
  case 0:
    variance = val;
    break;
  default:
    throw ndlexceptions::Error("Requested parameter doesn't exist.");
  }
}
double CLinKern::getParam(unsigned int paramNo) const
{
  BOUNDCHECK(paramNo==0);
  switch(paramNo)
  {
  case 0:
    return variance;
    break;
  default:
    throw ndlexceptions::Error("Requested parameter doesn't exist.");
  }
}
void CLinKern::getGradX(CMatrix& gX, const CMatrix& X, unsigned int row, const CMatrix& X2, bool addG) const
{
  DIMENSIONMATCH(gX.getRows() == X2.getRows());
  BOUNDCHECK(row < X.getRows());
  DIMENSIONMATCH(X.getCols()==X2.getCols());  

  for(unsigned int k=0; k<X2.getRows(); k++)
  {
    for(unsigned int j=0; j<X2.getCols(); j++)
    {
      double val = 0.0;
      if(addG)
	val = gX.getVal(k, j);
      val += variance*X2.getVal(k, j);
      gX.setVal(val, k, j);
    }
  }
}

void CLinKern::getDiagGradX(CMatrix& gX, const CMatrix& X, bool addG) const
{
  DIMENSIONMATCH(gX.dimensionsMatch(X));
  if(!addG)
  {
    gX.deepCopy(X);
    gX.scale(2.0*variance);
  }
  else
  {
    gX.axpy(X, 2.0*variance);
  }
}
double CLinKern::getWhite() const
{
  return 0.0;
}

double CLinKern::computeElement(const CMatrix& X1, unsigned int index1, 
				const CMatrix& X2, unsigned int index2) const
{
  return variance*X1.dotRowRow(index1, X2, index2);
}

void CLinKern::compute(CMatrix& K, const CMatrix& X) const
{
  DIMENSIONMATCH(K.rowsMatch(X));
  MATRIXPROPERTIES(K.isSquare());
  K.setSymmetric(true);
  K.syrk(X, variance, 0.0, "u", "n");
}

void CLinKern::compute(CMatrix& K, const CMatrix& X, const CMatrix& X2) const
{
  DIMENSIONMATCH(K.rowsMatch(X));
  DIMENSIONMATCH(K.getCols()==X2.getRows());
  K.gemm(X, X2, variance, 0.0, "n", "t");
}
void CLinKern::compute(CMatrix& K, const CMatrix& X, const CMatrix& X2, unsigned int row) const
{
  DIMENSIONMATCH(K.rowsMatch(X));
  DIMENSIONMATCH(K.getCols()==1);
  K.gemvRowRow(0, X, X2, row, variance, 0.0, "n");
}
double CLinKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad) const
{
  BOUNDCHECK(index==0);
  DIMENSIONMATCH(X.getRows()==covGrad.getRows());
  DIMENSIONMATCH(X2.getRows()==covGrad.getCols());
  double dot;
  double g1=0.0;
  for(unsigned int i=0; i<X.getRows(); i++)
    for(unsigned int j=0; j<X2.getRows(); j++)
    {
      dot = X.dotRowRow(i, X2, j);
      g1 += dot*covGrad.getVal(i, j);
    }
  return g1;
}
double CLinKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& covGrad) const
{
  BOUNDCHECK(index==0);
  DIMENSIONMATCH(X.rowsMatch(covGrad));
  MATRIXPROPERTIES(covGrad.isSquare());
  double dot;
  double g1=0.0;
  for(unsigned int i=0; i<X.getRows(); i++)
    for(unsigned int j=0; j<X.getRows(); j++)
    {
      dot = X.dotRowRow(i, X, j);
      g1 += dot*covGrad.getVal(i, j);
    }
  return g1;
}
  

CMlpKern::CMlpKern() : CKern()
{
  _init();
}
CMlpKern::CMlpKern(unsigned int inDim) : CKern(inDim)
{
  _init();
  setInputDim(inDim);
}
CMlpKern::CMlpKern(const CMatrix& X) : CKern(X)
{
  _init();
  setInputDim(X.getCols());
}  
CMlpKern::CMlpKern(const CMlpKern& kern) : CKern(kern)
{
  _init();
  setInputDim(kern.getInputDim());
  variance = kern.variance;
  weightVariance = kern.weightVariance;
  biasVariance = kern.biasVariance;
}
// Class destructor
CMlpKern::~CMlpKern()
{
}
double CMlpKern::getVariance() const
{
  return variance;
}
void CMlpKern::_init()
{
  nParams = 3;
  setType("mlp");
  setName("MLP");
  setParamName("weightVariance", 0);
  addTransform(CTransform::defaultPositive(), 0);
  setParamName("biasVariance", 1);
  addTransform(CTransform::defaultPositive(), 1);
  setParamName("variance", 2);
  addTransform(CTransform::defaultPositive(), 2);
  setStationary(false);
}
void CMlpKern::setInitParam()
{
  weightVariance = 10.0;
  biasVariance = 10.0;
  variance = 1.0;  
}

inline double CMlpKern::diagComputeElement(const CMatrix& X, unsigned int index) const
{
  double numer=weightVariance*X.norm2Row(index)+biasVariance;  
  double denom = numer+1.0;
  return variance*asin(numer/denom);
}
// Parameters are kernel parameters
void CMlpKern::setParam(double val, unsigned int paramNo)
{
  BOUNDCHECK(paramNo < nParams);
  switch(paramNo)
  {
  case 0:
    weightVariance = val;
    break;
  case 1:
    biasVariance = val;
    break;
  case 2:
    variance = val;
    break;
  default:
    throw ndlexceptions::Error("Requested parameter doesn't exist.");
  }
}
double CMlpKern::getParam(unsigned int paramNo) const
{
  BOUNDCHECK(paramNo < nParams);
  switch(paramNo)
  {
  case 0:
    return weightVariance;
    break;
  case 1:
    return biasVariance;
    break;
  case 2:
    return variance;
    break;
  default:
    throw ndlexceptions::Error("Requested parameter doesn't exist.");
  }
}

void CMlpKern::getGradX(CMatrix& gX, const CMatrix& X, unsigned int row, const CMatrix& X2, bool addG) const
{
  DIMENSIONMATCH(gX.getRows() == X2.getRows());
  BOUNDCHECK(row < X.getRows());
  DIMENSIONMATCH(X.getCols()==X2.getCols());
  double denomPart1 = weightVariance*X.norm2Row(row) + biasVariance + 1.0;
  for(unsigned int k=0; k<X2.getRows(); k++)
  {
    double numer=weightVariance*X.dotRowRow(row, X2, k) + biasVariance;
    double denomPart2=weightVariance*X2.norm2Row(k) + biasVariance + 1.0;
    double denom=sqrt(denomPart1*denomPart2);
    double arg = numer/denom;
    double kval=variance*weightVariance/sqrt(1-arg*arg);
    
    for(unsigned int j=0; j<X2.getCols(); j++)
    {
      double val = 0.0;
      if(addG)
	val = gX.getVal(k, j);
      val += (X2.getVal(k, j)/denom-X.getVal(row, j)*denomPart2*numer/(denom*denom*denom))*kval;
      gX.setVal(val, k, j);
    }
  } 
}
void CMlpKern::getDiagGradX(CMatrix& gX, const CMatrix& X, bool addG) const
{
  DIMENSIONMATCH(gX.dimensionsMatch(X));
  for(unsigned int i=0; i<X.getRows(); i++)
  {
    double numer=weightVariance*X.norm2Row(i) + biasVariance;
    double denom=numer+1.0;
    double arg=numer/denom;
    double kval=variance*weightVariance/sqrt(1-arg*arg);
    for(unsigned int j=0; j<X.getCols(); j++)
    {
      double val = 0.0;
      if(addG)
	val = gX.getVal(i, j);
      val += 2*(1/denom - numer/(denom*denom))*kval*X.getVal(i, j);
      gX.setVal(val, i, j);
    }
  }
}
double CMlpKern::getWhite() const
{
  return 0.0;
}

double CMlpKern::computeElement(const CMatrix& X1, unsigned int index1, 
				const CMatrix& X2, unsigned int index2) const
{
  double numer= weightVariance*X1.dotRowRow(index1, X2, index2) + biasVariance;
  double denom1=weightVariance*X1.norm2Row(index1)+biasVariance+1.0;  
  double denom2=weightVariance*X2.norm2Row(index2)+biasVariance+1.0;  
  return variance*asin(numer/sqrt(denom1*denom2));
}
void CMlpKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==nParams);
  double g1=0.0;
  double g2=0.0;
  double g3=0.0;
  innerProdi.resize(1, X.getRows());
  innerProdj.resize(1, X2.getRows());
  for(unsigned int j=0; j<X2.getRows(); j++)
    innerProdj.setVal(X2.norm2Row(j), j);
  // do off diagonal gradients first.
  for(unsigned int i=0; i<X.getRows(); i++)
  {
    innerProdi.setVal(X.norm2Row(i), i);
    for(unsigned int j=0; j<X2.getRows(); j++)
    {	  
      double crossProd=X.dotRowRow(i, X2, j);
      double numer=weightVariance*crossProd + biasVariance;
      double denomi=weightVariance*innerProdi.getVal(i)+biasVariance+1.0;  
      double denomj=weightVariance*innerProdj.getVal(j)+biasVariance+1.0;  
      double denom = sqrt(denomi*denomj);
      double arg = numer/denom;
      double baseCovGrad = variance/sqrt(1-arg*arg)*covGrad.getVal(i, j);
      double denom3=denom*denom*denom;
      g1+=baseCovGrad*(crossProd/denom-0.5*numer/denom3*(denomi*innerProdj.getVal(j) + innerProdi.getVal(i)*denomj));
      g2+=baseCovGrad*(1.0/denom-0.5*numer/denom3*((innerProdi.getVal(i)+innerProdj.getVal(j))*weightVariance +2.0*biasVariance+2.0));
      g3+=asin(arg)*covGrad.getVal(i, j);
    }
  }
  g.setVal(g1, 0);
  g.setVal(g2, 1);
  g.setVal(g3, 2);
  if(regularise)
    addPriorGrad(g);
  
}
void CMlpKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& covGrad, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==nParams);
  double g1=0.0;
  double g2=0.0;
  double g3=0.0;
  innerProdi.resize(1, X.getRows());

  // do off diagonal gradients first.
  for(unsigned int i=0; i<X.getRows(); i++)
  {
    innerProdi.setVal(X.norm2Row(i), i);
    for(unsigned int j=0; j<i; j++)
    {	  
      double crossProd=X.dotRowRow(i, X, j);
      double numer=weightVariance*crossProd + biasVariance;
      double denomi=weightVariance*innerProdi.getVal(i)+biasVariance+1.0;  
      double denomj=weightVariance*innerProdi.getVal(j)+biasVariance+1.0;  
      double denom = sqrt(denomi*denomj);
      double arg = numer/denom;
      double baseCovGrad = variance/sqrt(1-arg*arg)*covGrad.getVal(i, j);
      double denom3=denom*denom*denom;
      g1+=baseCovGrad*(crossProd/denom-0.5*numer/denom3*(denomi*innerProdi.getVal(j) + innerProdi.getVal(i)*denomj));
      g2+=baseCovGrad*(1.0/denom-0.5*numer/denom3*((innerProdi.getVal(i)+innerProdi.getVal(j))*weightVariance +2.0*biasVariance+2.0));
      g3+=asin(arg)*covGrad.getVal(i, j);
    }
  }
  // double the result due to symmetry.
  g1*=2;
  g2*=2;
  g3*=2;
  // add effect of diagonals.
  for(unsigned int i=0; i<X.getRows(); i++)
  {
    double numer = weightVariance*innerProdi.getVal(i)+biasVariance;
    double denom = numer+1.0;
    double denom3=denom*denom*denom;
    double arg = numer/denom;
    double baseCovGrad = variance/sqrt(1-arg*arg)*covGrad.getVal(i, i);
    g1+=baseCovGrad*(innerProdi.getVal(i)/denom-0.5*numer/denom3*(2*denom*innerProdi.getVal(i)));
      
    g2+=baseCovGrad*(1.0/denom-0.5*numer/denom3*(2.0*weightVariance*innerProdi.getVal(i)+2.0*biasVariance+2.0));
    g3+=asin(arg)*covGrad.getVal(i, i);

  }
  g.setVal(g1, 0);
  g.setVal(g2, 1);
  g.setVal(g3, 2);
  if(regularise)
    addPriorGrad(g);

}
double CMlpKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad) const
{
  
  BOUNDCHECK(index<nParams);
  throw ndlexceptions::NotImplementedError( "Error getGradParam is not currently implemented for CMlpKern");
  
}
double CMlpKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& covGrad) const
{
  
  BOUNDCHECK(index<nParams);
  throw ndlexceptions::NotImplementedError( "Error getGradParam is not currently implemented for CMlpKern");
  
}

CPolyKern::CPolyKern() : CKern()
{
  _init();
}
CPolyKern::CPolyKern(unsigned int inDim) : CKern(inDim)
{
  _init();
  setInputDim(inDim);
}
CPolyKern::CPolyKern(const CMatrix& X) : CKern(X)
{
  _init();
  setInputDim(X.getCols());
}  
CPolyKern::CPolyKern(const CPolyKern& kern) : CKern(kern)
{
  _init();
  setInputDim(kern.getInputDim());
  variance = kern.variance;
  weightVariance = kern.weightVariance;
  biasVariance = kern.biasVariance;
  degree = kern.degree;
}
// Class destructor
CPolyKern::~CPolyKern()
{
}
void CPolyKern::writeParamsToStream(ostream& out) const
{
  writeToStream(out, "baseType", getBaseType());
  writeToStream(out, "type", getType());
  writeToStream(out, "inputDim", getInputDim());
  writeToStream(out, "numParams", getNumParams());
  double deg = getDegree();
  if((deg - (int)deg)==0)
    writeToStream(out, "degree", (int)deg);
  else
    writeToStream(out, "degree", deg);
  CMatrix par(1, getNumParams());
  getParams(par);
  par.toStream(out);
  writeToStream(out, "numPriors", getNumPriors());
  writePriorsToStream(out);
}
void CPolyKern::readParamsFromStream(istream& in)
{
//   string tbaseType = getBaseTypeStream(in);
//   if(tbaseType != getBaseType())
//     throw ndlexceptions::StreamFormatError("baseType", "Error mismatch between saved base type, " + tbaseType + ", and Class base type, " + getType() + ".");
//   string ttype = getTypeStream(in);
//   if(ttype != getType())
//     throw ndlexceptions::StreamFormatError("type", "Error mismatch between saved type, " + ttype + ", and Class type, " + getType() + ".");
  setInputDim(readIntFromStream(in, "inputDim"));
  unsigned int numParams=readIntFromStream(in, "numParams");
  setDegree(readDoubleFromStream(in, "degree"));

  CMatrix par(1, numParams);
  par.fromStream(in);
  if(numParams==getNumParams())
    setParams(par);
  else
    throw ndlexceptions::StreamFormatError("numParams", "Listed number of parameters does not match computed number of parameters.");
  unsigned int numPriors = readIntFromStream(in, "numPriors");
  readPriorsFromStream(in, numPriors);
}
double CPolyKern::getVariance() const
{
  return variance;
}
void CPolyKern::_init()
{
  nParams = 3;
  setType("poly");
  setName("Polynomial");
  setParamName("weightVariance", 0);
  addTransform(CTransform::defaultPositive(), 0);
  setParamName("biasVariance", 1);
  addTransform(CTransform::defaultPositive(), 1);
  setParamName("variance", 2);
  addTransform(CTransform::defaultPositive(), 2);
  setStationary(false);
}
void CPolyKern::setInitParam()
{
  weightVariance = 1.0;
  biasVariance = 1.0;
  variance = 1.0;
  degree = 2.0;
}

inline double CPolyKern::diagComputeElement(const CMatrix& X, unsigned int index) const
{
  double arg=weightVariance*X.norm2Row(index)+biasVariance;  
  return variance*pow(arg, degree);
}
// Parameters are kernel parameters
void CPolyKern::setParam(double val, unsigned int paramNo)
{
  BOUNDCHECK(paramNo < nParams);
  switch(paramNo)
  {
  case 0:
    weightVariance = val;
    break;
  case 1:
    biasVariance = val;
    break;
  case 2:
    variance = val;
    break;
  default:
    throw ndlexceptions::Error("Requested parameter doesn't exist.");
  }
}
double CPolyKern::getParam(unsigned int paramNo) const
{
  BOUNDCHECK(paramNo < nParams);
  switch(paramNo)
  {
  case 0:
    return weightVariance;
    break;
  case 1:
    return biasVariance;
    break;
  case 2:
    return variance;
    break;
  default:
    throw ndlexceptions::Error("Requested parameter doesn't exist.");
  }
}

void CPolyKern::getGradX(CMatrix& gX, const CMatrix& X, unsigned int row, const CMatrix& X2, bool addG) const
{
  DIMENSIONMATCH(gX.getRows() == X2.getRows());
  BOUNDCHECK(row < X.getRows());
  DIMENSIONMATCH(X.getCols()==X2.getCols());
  for(unsigned int k=0; k<X2.getRows(); k++)
  {
    double arg=weightVariance*X.dotRowRow(row, X2, k) + biasVariance;
    double kval=degree*variance*weightVariance*pow(arg, degree-1);
    for(unsigned int j=0; j<X2.getCols(); j++)
    {
      double val = 0.0;
      if(addG)
	val = gX.getVal(k, j);
      val += kval*X2.getVal(k, j);
      gX.setVal(val, k, j);
    }
  }
}
void CPolyKern::getDiagGradX(CMatrix& gX, const CMatrix& X, bool addG) const
{
  DIMENSIONMATCH(gX.dimensionsMatch(X));
  for(unsigned int i=0; i<X.getRows(); i++)
  {
    double arg=weightVariance*X.norm2Row(i) + biasVariance;
    double kval=degree*variance*weightVariance*pow(arg, degree-1);
    for(unsigned int j=0; j<X.getCols(); j++)
    {
      double val = 0.0;
      if(addG)
	val = gX.getVal(i, j);
      val += 2*kval*X.getVal(i, j);
      gX.setVal(val, i, j);
    }
  }
}
double CPolyKern::getWhite() const
{
  return 0.0;
}

double CPolyKern::computeElement(const CMatrix& X1, unsigned int index1, 
				 const CMatrix& X2, unsigned int index2) const
{
  double arg=weightVariance*X1.dotRowRow(index1, X2, index2) + biasVariance;
  return variance*pow(arg, degree);
}
void CPolyKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==nParams);
  double g1=0.0;
  double g2=0.0;
  double g3=0.0;
  innerProdi.resize(1, X.getRows());
  // do off diagonal gradients first.
  for(unsigned int i=0; i<X.getRows(); i++) {
    innerProdi.setVal(X.norm2Row(i),i);
    for(unsigned int j=0; j<X2.getRows(); j++) {
      double crossProd=X.dotRowRow(i, X2, j);
      double arg=weightVariance*crossProd + biasVariance;
      double base = variance*degree*pow(arg, degree-1)*covGrad.getVal(i, j);	  
      g1+=crossProd*base;
      g2+=base;
      g3+=pow(arg, degree)*covGrad.getVal(i, j);
    }
  }
  g.setVal(g1, 0);
  g.setVal(g2, 1);
  g.setVal(g3, 2);
  if(regularise)
    addPriorGrad(g);

}
void CPolyKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& covGrad, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==nParams);
  double g1=0.0;
  double g2=0.0;
  double g3=0.0;
  innerProdi.resize(1, X.getRows());
  // do off diagonal gradients first.
  for(unsigned int i=0; i<X.getRows(); i++)
  {
    innerProdi.setVal(X.norm2Row(i),i);
    for(unsigned int j=0; j<i; j++)
    {
	  
      double crossProd=X.dotRowRow(i, X, j);
      double arg=weightVariance*crossProd + biasVariance;
      double base = variance*degree*pow(arg, degree-1)*covGrad.getVal(i, j);	  
      g1+=crossProd*base;
      g2+=base;
      g3+=pow(arg, degree)*covGrad.getVal(i, j);
    }
  }
  // double the result due to symmetry.
  g1*=2;
  g2*=2;
  g3*=2;
  // add effect of diagonals.
  for(unsigned int i=0; i<X.getRows(); i++)
  {
    double arg = weightVariance*innerProdi.getVal(i)+biasVariance;
    double base = variance*degree*pow(arg, degree-1)*covGrad.getVal(i, i);	  
    g1+=innerProdi.getVal(i)*base;
    g2+=base;
    g3+=pow(arg,degree)*covGrad.getVal(i, i);

  }
  g.setVal(g1, 0);
  g.setVal(g2, 1);
  g.setVal(g3, 2);
  if(regularise)
    addPriorGrad(g);

}
double CPolyKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad) const
{
  
  BOUNDCHECK(index<nParams);
  throw ndlexceptions::NotImplementedError( "Error getGradParam is not currently implemented for CPolyKern");
  
}
double CPolyKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& covGrad) const
{
  
  BOUNDCHECK(index<nParams);
  throw ndlexceptions::NotImplementedError( "Error getGradParam is not currently implemented for CPolyKern");
  
}
#ifdef _NDLMATLAB
void CPolyKern::addParamToMxArray(mxArray* matlabArray) const
{
  mxAddField(matlabArray, "degree");
  mxSetField(matlabArray, 0, "degree", convertMxArray(degree));
  CKern::addParamToMxArray(matlabArray);
}
void CPolyKern::extractParamFromMxArray(const mxArray* matlabArray)
{
  degree = mxArrayExtractIntField(matlabArray, "degree");
  CKern::extractParamFromMxArray(matlabArray);
}
#endif

// the Linear ARD kernel.
CLinardKern::CLinardKern() : CArdKern()
{
  _init();
}
CLinardKern::CLinardKern(unsigned int inDim) : CArdKern(inDim)
{
  _init();
  setInputDim(inDim);
}
CLinardKern::CLinardKern(const CMatrix& X) : CArdKern(X)
{
  _init();
  setInputDim(X.getCols());
}  
CLinardKern::CLinardKern(const CLinardKern& kern) : CArdKern(kern)
{
  _init();
  setInputDim(kern.getInputDim());
  variance = kern.variance;
  scales = kern.scales;
}
// Class destructor
CLinardKern::~CLinardKern()
{
}
double CLinardKern::getVariance() const
{
  return variance;
}
void CLinardKern::_init()
{
  nParams = 1;
  setType("linard");
  setName("linear ARD");
  setParamName("variance", 0);
  addTransform(CTransform::defaultPositive(), 0);
  setStationary(false);
}
void CLinardKern::setInitParam()
{ 
  nParams = 1+getInputDim();
  variance = 1.0;
  scales.resize(1, getInputDim());
  for(unsigned int i=1; i<getInputDim()+1; i++)
  {
    string name = "inputScale";
    setParamName(name, i);
  }
  scales.setVals(0.5);
  for(unsigned int i=1; i<getInputDim()+1; i++)
  {
    addTransform(CTransform::defaultZeroOne(), i);
  }
}

double CLinardKern::diagComputeElement(const CMatrix& X, unsigned int index1) const
{
  double val = 0.0;
  for(unsigned int i=0; i<getInputDim(); i++)
  {
    double x=X.getVal(index1, i);
    val += x*x*scales.getVal(i);
  }
  return val*variance;
}
// Parameters are kernel parameters
void CLinardKern::setParam(double val, unsigned int paramNo)
{
  BOUNDCHECK(paramNo<nParams);
  switch(paramNo)
  {
  case 0:
    variance = val;
    break;
  default:
    if(paramNo<getInputDim()+1)
      scales.setVal(val, paramNo-1);
    else
    {
      throw ndlexceptions::Error("Requested parameter doesn't exist.");
    }  
  }
}
double CLinardKern::getParam(unsigned int paramNo) const
{
  
  BOUNDCHECK(paramNo<nParams);
  switch(paramNo)
  {
  case 0:
    return variance;
    break;
  default:
    if(paramNo<getInputDim()+1)
      return scales.getVal(paramNo-1);
    else
    {
      throw ndlexceptions::Error("Requested parameter doesn't exist.");
    }  
  }
}
void CLinardKern::getGradX(CMatrix& gX, const CMatrix& X, unsigned int row, const CMatrix& X2, bool addG) const
{
  DIMENSIONMATCH(gX.getRows() == X2.getRows());
  BOUNDCHECK(row < X.getRows());
  DIMENSIONMATCH(X.getCols()==X2.getCols());
  for(unsigned int k=0; k<X2.getRows(); k++)
  {
    for(unsigned int j=0; j<X2.getCols(); j++)
    {
      double val = 0.0;
      if(addG)
	val = gX.getVal(k, j);
      val += variance*X2.getVal(k, j)*scales.getVal(j);
      gX.setVal(val, k, j);
    }
  }
}
void CLinardKern::getDiagGradX(CMatrix& gX, const CMatrix& X, bool addG) const
{
  DIMENSIONMATCH(gX.dimensionsMatch(X));
  if(!addG)
  {
    gX.deepCopy(X);
    double baseScale = 2.0*variance;
    for(unsigned int j=0; j<gX.getCols(); j++)
      gX.scaleCol(j, baseScale*scales.getVal(j));
  }
  else
  {
    double baseScale = 2.0*variance;
    for(unsigned int j=0; j<gX.getCols(); j++)
      gX.axpyColCol(j, X, j, baseScale*scales.getVal(j));
  }
}
double CLinardKern::getWhite() const
{
  return 0.0;
}

double CLinardKern::computeElement(const CMatrix& X1, unsigned int index1, 
				   const CMatrix& X2, unsigned int index2) const
{
  double val = 0.0;
  for(unsigned int i=0; i<getInputDim(); i++)
  {
    val += X1.getVal(index1, i)*X2.getVal(index2, i)*scales.getVal(i);
  }
  return val*variance;
}
void CLinardKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==nParams);
  double g1=0.0;
  for(unsigned int i=0; i<X.getRows(); i++) {
    for(unsigned int j=0; j<X2.getRows(); j++) {
      double val=0.0;
      for(unsigned int k=0; k<getInputDim(); k++)
	val += X.getVal(i, k)*X2.getVal(j, k)*scales.getVal(k);
      g1+=covGrad.getVal(i, j)*val;
    }
  }
  g.setVal(g1, 0);
  for(unsigned int k=0; k<getInputDim(); k++) {
    double g2=0.0;
    for(unsigned int i=0; i<X.getRows(); i++) {
      for(unsigned int j=0; j<X2.getRows(); j++) {
	g2+=X.getVal(i, k)*X2.getVal(j, k)*covGrad.getVal(i, j);
      }
    }
    g2*=variance;
    g.setVal(g2, k+1);
  }
  if(regularise)
    addPriorGrad(g);

}
void CLinardKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& covGrad, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==nParams);
  double g1=0.0;
  for(unsigned int i=0; i<X.getRows(); i++)
  {
    for(unsigned int j=0; j<i; j++)
    {
      double val=0.0;
      for(unsigned int k=0; k<getInputDim(); k++)
	val += X.getVal(i, k)*X.getVal(j, k)*scales.getVal(k);
      g1+=covGrad.getVal(i, j)*val;
    }
  }
  g1*=2.0;
  for(unsigned int i=0; i<X.getRows(); i++)
  {
    double val=0.0;
    for(unsigned int k=0; k<getInputDim(); k++)
      val += X.getVal(i, k)*X.getVal(i, k)*scales.getVal(k);
    g1+=covGrad.getVal(i, i)*val;
  }
  g.setVal(g1, 0);
  for(unsigned int k=0; k<getInputDim(); k++)
  {
    double g2=0.0;
    for(unsigned int i=0; i<X.getRows(); i++)
    {
      for(unsigned int j=0; j<i; j++)
      {
	g2+=X.getVal(i, k)*X.getVal(j, k)*covGrad.getVal(i, j);
      }
    }
    g2*=2;
    for(unsigned int i=0; i<X.getRows(); i++)
      g2+=X.getVal(i, k)*X.getVal(i, k)*covGrad.getVal(i, i);
    g2*=variance;
    g.setVal(g2, k+1);
  }
  if(regularise)
    addPriorGrad(g);

}
double CLinardKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad) const
{
  
  BOUNDCHECK(index<nParams);
  throw ndlexceptions::NotImplementedError( "Error getGradParam is not currently implemented for CLinardKern");
  
}
double CLinardKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& covGrad) const
{
  
  BOUNDCHECK(index<nParams);
  throw ndlexceptions::NotImplementedError( "Error getGradParam is not currently implemented for CLinardKern");
  
}
// the RBF ARD kernel.
CRbfardKern::CRbfardKern() : CArdKern()
{
  _init();
}
CRbfardKern::CRbfardKern(unsigned int inDim) : CArdKern(inDim)
{
  _init();
  setInputDim(inDim);
}
CRbfardKern::CRbfardKern(const CMatrix& X) : CArdKern(X)
{
  _init();
  setInputDim(X.getCols());
}  
CRbfardKern::CRbfardKern(const CRbfardKern& kern) : CArdKern(kern)
{
  _init();
  setInputDim(kern.getInputDim());
  variance = kern.variance;
  inverseWidth = kern.inverseWidth;
  scales = kern.scales;
}
// Class destructor
CRbfardKern::~CRbfardKern()
{
}
double CRbfardKern::getVariance() const
{
  return variance;
}
void CRbfardKern::_init()
{
  nParams = 2;
  setType("rbfard");
  setName("RBF ARD");
  setParamName("inverseWidth", 0);
  addTransform(CTransform::defaultPositive(), 0);
  setParamName("variance", 1);
  addTransform(CTransform::defaultPositive(), 1);
  setStationary(true);
}
void CRbfardKern::setInitParam()
{
  nParams = 2+getInputDim();
  inverseWidth=1.0;
  variance = 1.0;

  // input scales.
  scales.resize(1, getInputDim());
  gscales.resize(1, getInputDim());
  scales.setVals(0.5);
  for(unsigned int i=2; i<getInputDim()+2; i++)
  {
    string name = "inputScale";
    setParamName(name, i);
  }
  for(unsigned int i=2; i<getInputDim()+2; i++)
  {
    addTransform(CTransform::defaultZeroOne(), i);
  }

}

double CRbfardKern::diagComputeElement(const CMatrix& X, unsigned int index1) const
{
  return variance;
}
// Parameters are kernel parameters
void CRbfardKern::setParam(double val, unsigned int paramNo)
{
  
  BOUNDCHECK(paramNo<nParams);
  switch(paramNo)
  {
  case 0:
    inverseWidth=val;
    break;
  case 1:
    variance=val;
    break;
  default:
    if(paramNo<nParams)
      scales.setVal(val, paramNo-2);
    else
    {
      throw ndlexceptions::Error("Requested parameter doesn't exist.");
    }    
  }
}
double CRbfardKern::getParam(unsigned int paramNo) const
{
  
  BOUNDCHECK(paramNo<nParams);
  switch(paramNo)
  {
  case 0:
    return inverseWidth;
    break;
  case 1:
    return variance;
    break;
  default:
    if(paramNo<nParams)
      return scales.getVal(paramNo-2);
    else
    {
      throw ndlexceptions::Error("Requested parameter doesn't exist.");
    }  
  }
}
void CRbfardKern::getGradX(CMatrix& gX, const CMatrix& X, unsigned int row, const CMatrix& X2, bool addG) const
{
  DIMENSIONMATCH(gX.getRows() == X2.getRows());
  BOUNDCHECK(row < X.getRows());
  DIMENSIONMATCH(X.getCols()==X2.getCols());
  double wi2 = 0.5*inverseWidth;
  double pf = variance*inverseWidth;
  for(unsigned int k=0; k<X2.getRows(); k++)
  {
    double n2=0.0;
    for(unsigned int j=0; j<getInputDim(); j++)
    {
      double x = X.getVal(row, j);
      x = x-X2.getVal(k, j);
      n2+=x*scales.getVal(j)*x;
    }
    for(unsigned int j=0; j<X2.getCols(); j++)
    {
      double val = 0.0;
      if(addG)
	val = gX.getVal(k, j);
      val += pf*(X2.getVal(k, j)-X.getVal(row, j))*exp(-n2*wi2)*scales.getVal(j);
      gX.setVal(val, k, j);	      
    }
  }
}
void CRbfardKern::getDiagGradX(CMatrix& gX, const CMatrix& X, bool addG) const
{
  DIMENSIONMATCH(gX.dimensionsMatch(X));
  if(!addG)
    gX.zeros();
}
double CRbfardKern::getWhite() const
{
  return 0.0;
}

double CRbfardKern::computeElement(const CMatrix& X1, unsigned int index1, 
				   const CMatrix& X2, unsigned int index2) const
{
  double val = 0.0;
  for(unsigned int i=0; i<getInputDim(); i++)
  {
    double x = X1.getVal(index1, i);
    x = x-X2.getVal(index2, i);
    val+=x*scales.getVal(i)*x;
  }
  return variance*exp(-val*inverseWidth*0.5);
}

void CRbfardKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==nParams);
  double g1=0.0;
  double g2=0.0;
  gscales.zeros();
  double halfInverseWidth = 0.5*inverseWidth;
  for(unsigned int i=0; i<X.getRows(); i++)
  {
    for(unsigned int j=0; j<X2.getRows(); j++)
    {
      double val = 0.0;
      for(unsigned int k=0; k<getInputDim(); k++)
      {
	double x = X.getVal(i, k);
	x-=X2.getVal(j, k);
	val+=x*scales.getVal(k)*x;
      }
      double kCovGrad = exp(-halfInverseWidth*val)*covGrad.getVal(i, j);
      g1-=0.5*val*kCovGrad*variance;
      g2+=kCovGrad;
      for(unsigned int k=0; k<getInputDim(); k++)
      {
	double g3=gscales.getVal(k);
	double xi=X.getVal(i, k);
	double xj=X2.getVal(j, k);
	g3+=inverseWidth*kCovGrad*(xi*xj-.5*xi*xi-.5*xj*xj)*variance;
	gscales.setVal(g3, k);
      }
    }
  }
  g.setVal(g1, 0);
  g.setVal(g2, 1);
  for(unsigned int k=0; k<getInputDim(); k++)
    g.setVal(gscales.getVal(k), k+2);
  if(regularise)
    addPriorGrad(g);

}

void CRbfardKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& covGrad, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==nParams);
  double g1=0.0;
  double g2=0.0;
  gscales.zeros();
  double halfInverseWidth = 0.5*inverseWidth;
  for(unsigned int i=0; i<X.getRows(); i++)
  {
    for(unsigned int j=0; j<i; j++)
    {
      double val = 0.0;
      for(unsigned int k=0; k<getInputDim(); k++)
      {
	double x = X.getVal(i, k);
	x-=X.getVal(j, k);
	val+=x*scales.getVal(k)*x;
      }
      double kCovGrad = exp(-halfInverseWidth*val)*covGrad.getVal(i, j);
      g1-=0.5*val*kCovGrad*variance;
      g2+=kCovGrad;
      for(unsigned int k=0; k<getInputDim(); k++)
      {
	double g3=gscales.getVal(k);
	double xi=X.getVal(i, k);
	double xj=X.getVal(j, k);
	g3+=inverseWidth*kCovGrad*(xi*xj-.5*xi*xi-.5*xj*xj)*variance;
	gscales.setVal(g3, k);
      }
    }
  }
  g1*=2.0;
  g2*=2.0;
  for(unsigned int i=0; i<X.getRows(); i++)
    g2+=covGrad.getVal(i, i);
  gscales.scale(2.0);
  g.setVal(g1, 0);
  g.setVal(g2, 1);
  for(unsigned int k=0; k<getInputDim(); k++)
    g.setVal(gscales.getVal(k), k+2);
  if(regularise)
    addPriorGrad(g);

}
double CRbfardKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad) const
{
  
  BOUNDCHECK(index<nParams);
  throw ndlexceptions::NotImplementedError( "Error getGradParam is not currently implemented for CRbfardKern");
  
}
double CRbfardKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& covGrad) const
{
  
  BOUNDCHECK(index<nParams);
  throw ndlexceptions::NotImplementedError( "Error getGradParam is not currently implemented for CRbfardKern");
  
}

// the MLP ARD kernel.
CMlpardKern::CMlpardKern() : CArdKern()
{
}
CMlpardKern::CMlpardKern(unsigned int inDim) : CArdKern(inDim)
{
  _init();
  setInputDim(inDim);
}
CMlpardKern::CMlpardKern(const CMatrix& X) : CArdKern(X) 
{
  _init();
  setInputDim(X.getCols());
}  
CMlpardKern::CMlpardKern(const CMlpardKern& kern) : CArdKern(kern)
{
  _init();
  setInputDim(kern.getInputDim());
  variance = kern.variance;
  weightVariance=kern.weightVariance;
  biasVariance=kern.biasVariance;
  scales = kern.scales;
}
// Class destructor
CMlpardKern::~CMlpardKern()
{
}
double CMlpardKern::getVariance() const
{
  return variance;
}
void CMlpardKern::_init()
{
  nParams = 3;
  setType("mlpard");
  setName("MLP ARD");
  setParamName("weightVariance", 0);
  addTransform(CTransform::defaultPositive(), 0);
  setParamName("biasVariance", 1);
  addTransform(CTransform::defaultPositive(), 1);
  setParamName("variance", 2);
  addTransform(CTransform::defaultPositive(), 2);
  setStationary(false);

}
void CMlpardKern::setInitParam()
{
  nParams = 3+getInputDim();
  weightVariance=10.0;
  biasVariance=10.0;
  variance = 1.0;
  scales.resize(1, getInputDim());
  gscales.resize(1, getInputDim());
  scales.setVals(0.5);
  for(unsigned int i=3; i<getInputDim()+3; i++)
  {
    string name = "inputScale";
    setParamName(name, i);
  }
  for(unsigned int i=3; i<getInputDim()+3; i++)
  {
    addTransform(CTransform::defaultZeroOne(), i);
  }
}

double CMlpardKern::diagComputeElement(const CMatrix& X, unsigned int index) const
{
  double val=0.0;
  for(unsigned int i=0; i<getInputDim(); i++)
  {
    double x = X.getVal(index, i);
    val+=x*x*scales.getVal(i);
  }
  double numer=weightVariance*val+biasVariance;  
  double denom = numer+1.0;
  return variance*asin(numer/denom);
}
// Parameters are kernel parameters
void CMlpardKern::setParam(double val, unsigned int paramNo)
{
  
  BOUNDCHECK(paramNo<nParams);
  switch(paramNo)
  {
  case 0:
    weightVariance=val;
    break;
  case 1:
    biasVariance=val;
    break;
  case 2:
    variance=val;
    break;
  default:
    if(paramNo<nParams)
      scales.setVal(val, paramNo-3);
    else
    {
      throw ndlexceptions::Error("Requested parameter doesn't exist.");
    }  
  }
}
double CMlpardKern::getParam(unsigned int paramNo) const
{
  
  BOUNDCHECK(paramNo<nParams);
  switch(paramNo)
  {
  case 0:
    return weightVariance;
    break;
  case 1:
    return biasVariance;
    break;
  case 2:
    return variance;
    break;
  default:
    if(paramNo<nParams)
      return scales.getVal(paramNo-3);
    else
    {
      throw ndlexceptions::Error("Requested parameter doesn't exist.");
    }  
  }
}

void CMlpardKern::getGradX(CMatrix& gX, const CMatrix& X, unsigned int row, const CMatrix& X2, bool addG) const
{
  DIMENSIONMATCH(gX.getRows() == X2.getRows());
  BOUNDCHECK(row < X.getRows());
  DIMENSIONMATCH(X.getCols()==X2.getCols());
  double val=0.0;
  for(unsigned int j=0; j<getInputDim(); j++)
  {
    double x = X.getVal(row, j);
    val+=x*x*scales.getVal(j);
  }
  double denomPart1=weightVariance*val+biasVariance+1.0;  
  for(unsigned int k=0; k<X2.getRows(); k++)
  {
    double val=0.0;
    double val2=0.0;
    for(unsigned int j=0; j<getInputDim(); j++)
    {
      double xk = X2.getVal(k, j);
      double xi = X.getVal(row, j);
      double scalesX = xk*scales.getVal(j);
      val+=xk*scalesX;
      val2+=xi*scalesX;
    }
    double numer=weightVariance*val2 + biasVariance;
    double denomPart2=weightVariance*val+biasVariance+1.0;  
    double denom=sqrt(denomPart1*denomPart2);
    double arg = numer/denom;
    double kval=variance*weightVariance/sqrt(1-arg*arg);
    
    for(unsigned int j=0; j<X2.getCols(); j++)
    {
      double val = 0.0;
      if(addG)
	val = gX.getVal(k, j);
      val += (X2.getVal(k, j)/denom-X.getVal(row, j)*denomPart2*numer/(denom*denom*denom))*kval*scales.getVal(j);
      gX.setVal(val, k, j);
    }
  }
}
void CMlpardKern::getDiagGradX(CMatrix& gX, const CMatrix& X, bool addG) const
{
  DIMENSIONMATCH(gX.dimensionsMatch(X));
  for(unsigned int i=0; i<X.getRows(); i++)
  {
    double val=0.0;
    for(unsigned int j=0; j<getInputDim(); j++)
    {
      double x = X.getVal(i, j);
      val+=x*x*scales.getVal(j);
    }
    double numer=weightVariance*val+biasVariance;  
    double denom=numer+1.0;
    double arg=numer/denom;
    double kval=variance*weightVariance/sqrt(1-arg*arg);
    for(unsigned int j=0; j<X.getCols(); j++)
    {
      double val = 0.0;
      if(addG)
	val = gX.getVal(i, j);
      val += 2*(1/denom - numer/(denom*denom))*kval*X.getVal(i, j)*scales.getVal(j);
      gX.setVal(val, i, j);
    }
  }
}
double CMlpardKern::getWhite() const
{
  return 0.0;
}

double CMlpardKern::computeElement(const CMatrix& X1, unsigned int index1, 
				   const CMatrix& X2, unsigned int index2) const
{
  double valij=0.0;
  double valii=0.0;
  double valjj=0.0;
  for(unsigned int i=0; i<getInputDim(); i++)
  {
    double xi = X1.getVal(index1, i);
    double xj = X2.getVal(index2, i);
    valij+=xi*xj*scales.getVal(i);
    valii+=xi*xi*scales.getVal(i);
    valjj+=xj*xj*scales.getVal(i);
  }
  double numer=weightVariance*valij + biasVariance;
  double denom1=weightVariance*valii+biasVariance+1.0;  
  double denom2=weightVariance*valjj+biasVariance+1.0;  
  return variance*asin(numer/sqrt(denom1*denom2));
}

void CMlpardKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==nParams);
  double g1=0.0;
  double g2=0.0;
  double g3=0.0;
  innerProdi.resize(1, X.getRows());
  innerProdj.resize(1, X2.getRows());
  gscales.zeros();
  innerProdi.zeros();
  for(unsigned int i=0; i<X.getRows(); i++)
  {
    double val=0.0;
    for(unsigned int k=0; k<getInputDim(); k++)
    {
      double x = X.getVal(i, k);
      val+=x*x*scales.getVal(k);
    }
    innerProdi.setVal(val, i);
  }
  for(unsigned int j=0; j<X2.getRows(); j++)
  {
    double val=0.0;
    for(unsigned int k=0; k<getInputDim(); k++)
    {
      double x = X2.getVal(j, k);
      val+=x*x*scales.getVal(k);
    }
    innerProdj.setVal(val, j);
  }
  
  for(unsigned int i=0; i<X.getRows(); i++)
  {
    for(unsigned int j=0; j<X2.getRows(); j++)
    {
      double val=0.0;
      for(unsigned int k=0; k<getInputDim(); k++)
      {
	double xi = X.getVal(i, k);
	double xj = X2.getVal(j, k);
	val+=xi*xj*scales.getVal(k);
      }
      double numer=weightVariance*val + biasVariance;
      double denomi=weightVariance*innerProdi.getVal(i)+biasVariance+1.0;  
      double denomj=weightVariance*innerProdj.getVal(j)+biasVariance+1.0;  
      double denom = sqrt(denomi*denomj);
      double arg = numer/denom;
      double baseCovGrad = variance/sqrt(1-arg*arg)*covGrad.getVal(i, j);
      double denom3=denom*denom*denom;
      g1+=baseCovGrad*(val/denom-0.5*numer/denom3*(denomi*innerProdj.getVal(j) + innerProdi.getVal(i)*denomj));
      g2+=baseCovGrad*(1.0/denom-0.5*numer/denom3*((innerProdi.getVal(i)+innerProdj.getVal(j))*weightVariance+2.0*biasVariance+2.0));
      g3+=asin(arg)*covGrad.getVal(i, j);


      for(unsigned int k=0; k<getInputDim(); k++)
      {
	double g4=gscales.getVal(k);
	double xik=X.getVal(i, k);
	double xjk=X2.getVal(j, k);
	g4+=(xik*xjk/denom - 0.5*numer/denom3*(xik*xik*denomj + xjk*xjk*denomi))*baseCovGrad*weightVariance;
	gscales.setVal(g4, k);
      }
    }
  }
  g.setVal(g1, 0);
  g.setVal(g2, 1);
  g.setVal(g3, 2);
  for(unsigned int k=0; k<getInputDim(); k++)
    g.setVal(gscales.getVal(k), k+3);
  if(regularise)
    addPriorGrad(g);

}
void CMlpardKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& covGrad, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==nParams);
  double g1=0.0;
  double g2=0.0;
  double g3=0.0;
  innerProdi.resize(1, X.getRows());
  gscales.zeros();
  innerProdi.zeros();
  for(unsigned int i=0; i<X.getRows(); i++)
  {
    double val=0.0;
    for(unsigned int k=0; k<getInputDim(); k++)
    {
      double x = X.getVal(i, k);
      val+=x*x*scales.getVal(k);
    }
    innerProdi.setVal(val, i);
  }
  
  for(unsigned int i=0; i<X.getRows(); i++)
  {
    for(unsigned int j=0; j<i; j++)
    {
      double val=0.0;
      for(unsigned int k=0; k<getInputDim(); k++)
      {
	double xi = X.getVal(i, k);
	double xj = X.getVal(j, k);
	val+=xi*xj*scales.getVal(k);
      }
      double numer=weightVariance*val + biasVariance;
      double denomi=weightVariance*innerProdi.getVal(i)+biasVariance+1.0;  
      double denomj=weightVariance*innerProdi.getVal(j)+biasVariance+1.0;  
      double denom = sqrt(denomi*denomj);
      double arg = numer/denom;
      double baseCovGrad = variance/sqrt(1-arg*arg)*covGrad.getVal(i, j);
      double denom3=denom*denom*denom;
      g1+=baseCovGrad*(val/denom-0.5*numer/denom3*(denomi*innerProdi.getVal(j) + innerProdi.getVal(i)*denomj));
      g2+=baseCovGrad*(1.0/denom-0.5*numer/denom3*((innerProdi.getVal(i)+innerProdi.getVal(j))*weightVariance+2.0*biasVariance+2.0));
      g3+=asin(arg)*covGrad.getVal(i, j);


      for(unsigned int k=0; k<getInputDim(); k++)
      {
	double g4=gscales.getVal(k);
	double xik=X.getVal(i, k);
	double xjk=X.getVal(j, k);
	g4+=(xik*xjk/denom - 0.5*numer/denom3*(xik*xik*denomj + xjk*xjk*denomi))*baseCovGrad*weightVariance;
	gscales.setVal(g4, k);
      }
    }
  }
  // double due to symmetriy.
  g1*=2.0;
  g2*=2.0;
  g3*=2.0;
  gscales.scale(2.0);
  // add effect of diagonals.
  for(unsigned int i=0; i<X.getRows(); i++)
  {
    double numer = weightVariance*innerProdi.getVal(i)+biasVariance;
    double denom = numer+1.0;
    double denom3=denom*denom*denom;
    double arg = numer/denom;
    double baseCovGrad = variance/sqrt(1-arg*arg)*covGrad.getVal(i, i);
    g1+=baseCovGrad*(innerProdi.getVal(i)/denom-0.5*numer/denom3*(2*denom*innerProdi.getVal(i)));
      
    g2+=baseCovGrad*(1.0/denom-0.5*numer/denom3*(2.0*weightVariance*innerProdi.getVal(i)+2.0*biasVariance+2.0));
    g3+=asin(arg)*covGrad.getVal(i, i);

    for(unsigned int k=0; k<getInputDim(); k++)
    {
      double g4=gscales.getVal(k);
      double xik=X.getVal(i, k);
      g4+=(xik*xik/denom - 0.5*numer/denom3*(2*xik*xik*denom))*baseCovGrad*weightVariance;
      gscales.setVal(g4, k);
    }
  }
  g.setVal(g1, 0);
  g.setVal(g2, 1);
  g.setVal(g3, 2);
  for(unsigned int k=0; k<getInputDim(); k++)
    g.setVal(gscales.getVal(k), k+3);
  if(regularise)
    addPriorGrad(g);

}
double CMlpardKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& covGrad) const
{
  
  BOUNDCHECK(index<nParams);
  throw ndlexceptions::NotImplementedError( "Error getGradParam is not currently implemented for CMlpardKern");
  
}
double CMlpardKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad) const
{
  
  BOUNDCHECK(index<nParams);
  throw ndlexceptions::NotImplementedError( "Error getGradParam is not currently implemented for CMlpardKern");
  
}

// the POLY ARD kernel.
CPolyardKern::CPolyardKern() : CArdKern()
{
  _init();
}
CPolyardKern::CPolyardKern(unsigned int inDim) : CArdKern(inDim)
{
  _init();
  setInputDim(inDim);
  
}
CPolyardKern::CPolyardKern(const CMatrix& X) : CArdKern(X)
{
  _init();
  setInputDim(X.getCols());
}  
CPolyardKern::CPolyardKern(const CPolyardKern& kern) : CArdKern(kern)
{
  _init();
  setInputDim(kern.getInputDim());
  variance = kern.variance;
  weightVariance=kern.weightVariance;
  biasVariance=kern.biasVariance;
  scales = kern.scales;
}
// Class destructor
CPolyardKern::~CPolyardKern()
{
}
void CPolyardKern::writeParamsToStream(ostream& out) const
{
  writeToStream(out, "baseType", getBaseType());
  writeToStream(out, "type", getType());
  writeToStream(out, "inputDim", getInputDim());
  writeToStream(out, "numParams", getNumParams());
  double deg = getDegree();
  if((deg - (int)deg)==0)
    writeToStream(out, "degree", (int)deg);
  else
    writeToStream(out, "degree", deg);
  CMatrix par(1, getNumParams());
  getParams(par);
  par.toStream(out);
  writeToStream(out, "numPriors", getNumPriors());
  writePriorsToStream(out);



}
void CPolyardKern::readParamsFromStream(istream& in)
{
//   string tbaseType = getBaseTypeStream(in);
//   if(tbaseType != getBaseType())
//     throw ndlexceptions::StreamFormatError("baseType", "Error mismatch between saved base type, " + tbaseType + ", and Class base type, " + getType() + ".");
//   string ttype = getTypeStream(in);
//   if(ttype != getType())
//     throw ndlexceptions::StreamFormatError("type", "Error mismatch between saved type, " + ttype + ", and Class type, " + getType() + ".");
  setInputDim(readIntFromStream(in, "inputDim"));
  unsigned int numParams=readIntFromStream(in, "numParams");
  setDegree(readDoubleFromStream(in, "degree"));

  CMatrix par(1, numParams);
  par.fromStream(in);
  if(numParams==getNumParams())
    setParams(par);
  else
    throw ndlexceptions::StreamFormatError("numParams", "Listed number of parameters does not match computed number of parameters.");
  unsigned int numPriors = readIntFromStream(in, "numPriors");
  readPriorsFromStream(in, numPriors);

}
double CPolyardKern::getVariance() const
{
  return variance;
}
void CPolyardKern::_init()
{
  nParams = 3;
  setType("polyard");
  setName("Polynomial ARD");
  setParamName("weightVariance", 0);
  addTransform(CTransform::defaultPositive(), 0);
  setParamName("biasVariance", 1);
  addTransform(CTransform::defaultPositive(), 1);
  setParamName("variance", 2);
  addTransform(CTransform::defaultPositive(), 2);
  setStationary(false);
}
void CPolyardKern::setInitParam()
{
  nParams = 3+getInputDim();
  weightVariance=1.0;
  biasVariance=1.0;
  variance = 1.0;
  scales.resize(1, getInputDim());
  gscales.resize(1, getInputDim());
  scales.setVals(0.5);
  for(unsigned int i=3; i<getInputDim()+3; i++)
  {
    string name = "inputScale";
    setParamName(name, i);
  }
  for(unsigned int i=3; i<getInputDim()+3; i++)
  {
    addTransform(CTransform::defaultZeroOne(), i);
  }
  degree = 2.0;
}

double CPolyardKern::diagComputeElement(const CMatrix& X, unsigned int index) const
{
  double val=0.0;
  for(unsigned int i=0; i<getInputDim(); i++)
  {
    double x = X.getVal(index, i);
    val+=x*x*scales.getVal(i);
  }
  double arg=weightVariance*val+biasVariance;  
  return variance*pow(arg, degree);
}
// Parameters are kernel parameters
void CPolyardKern::setParam(double val, unsigned int paramNo)
{
  
  BOUNDCHECK(paramNo<nParams);
  switch(paramNo)
  {
  case 0:
    weightVariance=val;
    break;
  case 1:
    biasVariance=val;
    break;
  case 2:
    variance=val;
    break;
  default:
    if(paramNo<nParams)
      scales.setVal(val, paramNo-3);
    else
    {
      throw ndlexceptions::Error("Requested parameter doesn't exist.");
    }    
  }
}
double CPolyardKern::getParam(unsigned int paramNo) const
{
  
  BOUNDCHECK(paramNo<nParams);
  switch(paramNo)
  {
  case 0:
    return weightVariance;
    break;
  case 1:
    return biasVariance;
    break;
  case 2:
    return variance;
    break;
  default:
    if(paramNo<nParams)
      return scales.getVal(paramNo-3);
    else
    {
      throw ndlexceptions::Error("Requested parameter doesn't exist.");
    }  
  }
}
void CPolyardKern::getGradX(CMatrix& gX, const CMatrix& X, unsigned int row, const CMatrix& X2, bool addG) const
{
  DIMENSIONMATCH(gX.getRows() == X2.getRows());
  BOUNDCHECK(row < X.getRows());
  DIMENSIONMATCH(X.getCols()==X2.getCols());
  for(unsigned int k=0; k<X2.getRows(); k++)
  {
    double valik=0.0;
    for(unsigned int j=0; j<getInputDim(); j++)
    {
      valik+=X.getVal(row, j)*X2.getVal(k, j)*scales.getVal(j);
    }
    double arg=weightVariance*valik + biasVariance;
    double kval=degree*variance*weightVariance*pow(arg, degree-1);
    for(unsigned int j=0; j<X2.getCols(); j++)
    {
      double val = 0.0;
      if(addG)
	val = gX.getVal(k, j);
      val += kval*X2.getVal(k, j)*scales.getVal(j);
      gX.setVal(val, k, j);
    }
  }
}
void CPolyardKern::getDiagGradX(CMatrix& gX, const CMatrix& X, bool addG) const
{
  DIMENSIONMATCH(gX.dimensionsMatch(X));
  for(unsigned int i=0; i<X.getRows(); i++)
  {
    double val=0.0;
    for(unsigned int j=0; j<getInputDim(); j++)
    {
      double x = X.getVal(i, j);
      val+=x*x*scales.getVal(j);
    }
    double arg=weightVariance*val + biasVariance;
    double kval=degree*variance*weightVariance*pow(arg, degree-1);
    for(unsigned int j=0; j<X.getCols(); j++)
    {
      double val = 0.0;
      if(addG)
	val = gX.getVal(i, j);
      val += 2*kval*X.getVal(i, j)*scales.getVal(j);
      gX.setVal(val, i, j);
    }
  }
}
double CPolyardKern::getWhite() const
{
  return 0.0;
}

double CPolyardKern::computeElement(const CMatrix& X1, unsigned int index1, 
				    const CMatrix& X2, unsigned int index2) const
{
  double valij=0.0;
  for(unsigned int i=0; i<getInputDim(); i++)
  {
    valij+=X1.getVal(index1, i)*X2.getVal(index2, i)*scales.getVal(i);
  }
  double arg=weightVariance*valij + biasVariance;
  return variance*pow(arg, degree);
}

void CPolyardKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==nParams);
  double g1=0.0;
  double g2=0.0;
  double g3=0.0;
  innerProdi.resize(1, X.getRows());
  gscales.zeros();
  
  for(unsigned int i=0; i<X.getRows(); i++)
  {
    for(unsigned int j=0; j<X2.getRows(); j++)
    {
      double val=0.0;
      for(unsigned int k=0; k<getInputDim(); k++)
      {
	val+=X.getVal(i, k)*X2.getVal(j, k)*scales.getVal(k);
      }
      double arg=weightVariance*val + biasVariance;
      double baseCovGrad = variance*degree*pow(arg, degree-1)*covGrad.getVal(i, j);
      g1+=baseCovGrad*val;
      g2+=baseCovGrad;
      g3+=pow(arg, degree)*covGrad.getVal(i, j);


      for(unsigned int k=0; k<getInputDim(); k++)
      {
	double g4=gscales.getVal(k);
	g4+=X.getVal(i, k)*X2.getVal(j, k)*baseCovGrad*weightVariance;
	gscales.setVal(g4, k);
      }
    }
  }
  g.setVal(g1, 0);
  g.setVal(g2, 1);
  g.setVal(g3, 2);
  for(unsigned int k=0; k<getInputDim(); k++)
    g.setVal(gscales.getVal(k), k+3);
  if(regularise)
    addPriorGrad(g);

}
void CPolyardKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& covGrad, bool regularise) const
{
  DIMENSIONMATCH(g.getRows()==1);
  DIMENSIONMATCH(g.getCols()==nParams);
  double g1=0.0;
  double g2=0.0;
  double g3=0.0;
  innerProdi.resize(1, X.getRows());
  gscales.zeros();
  innerProdi.zeros();
  for(unsigned int i=0; i<X.getRows(); i++)
  {
    double val=0.0;
    for(unsigned int k=0; k<getInputDim(); k++)
    {
      double x = X.getVal(i, k);
      val+=x*x*scales.getVal(k);
    }
    innerProdi.setVal(val, i);
  }
  
  for(unsigned int i=0; i<X.getRows(); i++)
  {
    for(unsigned int j=0; j<i; j++)
    {
      double val=0.0;
      for(unsigned int k=0; k<getInputDim(); k++)
      {
	val+=X.getVal(i, k)*X.getVal(j, k)*scales.getVal(k);
      }
      double arg=weightVariance*val + biasVariance;
      double baseCovGrad = variance*degree*pow(arg, degree-1)*covGrad.getVal(i, j);
      g1+=baseCovGrad*val;
      g2+=baseCovGrad;
      g3+=pow(arg, degree)*covGrad.getVal(i, j);


      for(unsigned int k=0; k<getInputDim(); k++)
      {
	double g4=gscales.getVal(k);
	g4+=X.getVal(i, k)*X.getVal(j, k)*baseCovGrad*weightVariance;
	gscales.setVal(g4, k);
      }
    }
  }
  // double due to symmetry.
  g1*=2.0;
  g2*=2.0;
  g3*=2.0;
  gscales.scale(2.0);
  // add effect of diagonals.
  for(unsigned int i=0; i<X.getRows(); i++)
  {
    double arg = weightVariance*innerProdi.getVal(i)+biasVariance;
    double baseCovGrad = variance*degree*pow(arg, degree-1)*covGrad.getVal(i, i);
    g1+=baseCovGrad*innerProdi.getVal(i);
      
    g2+=baseCovGrad;
    g3+=pow(arg, degree)*covGrad.getVal(i, i);

    for(unsigned int k=0; k<getInputDim(); k++)
    {
      double g4=gscales.getVal(k);
      double xik=X.getVal(i, k);
      g4+=xik*xik*baseCovGrad*weightVariance;
      gscales.setVal(g4, k);
    }
  }
  g.setVal(g1, 0);
  g.setVal(g2, 1);
  g.setVal(g3, 2);
  for(unsigned int k=0; k<getInputDim(); k++)
    g.setVal(gscales.getVal(k), k+3);
  if(regularise)
    addPriorGrad(g);

}
double CPolyardKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& covGrad) const
{
  
  BOUNDCHECK(index<nParams);
  throw ndlexceptions::NotImplementedError( "Error getGradParam is not currently implemented for CPolyardKern");
  
}
double CPolyardKern::getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& covGrad) const
{
  
  BOUNDCHECK(index<nParams);
  throw ndlexceptions::NotImplementedError( "Error getGradParam is not currently implemented for CPolyardKern");
  
}

// Functions that operate on CKern.
ostream& operator<<(ostream& out, const CKern& kern)
{
  out <<  kern.display(out);
  return out;
}
void writeKernToStream(const CKern& kern, ostream& out)
{
  kern.toStream(out);
}
CKern* readKernFromStream(istream& in)
{
  double ver = CStreamInterface::readVersionFromStream(in); 
  string tbaseType = CStreamInterface::getBaseTypeStream(in);
  if(tbaseType != "kern")
    throw ndlexceptions::StreamFormatError("baseType", "Error mismatch between saved base type, " + tbaseType + ", and Class base type, kern.");
  CKern* pkern;

  string type = CStreamInterface::getTypeStream(in);
  if(type=="white")
    pkern = new CWhiteKern();
  else if(type=="whitefixed")
    pkern = new CWhitefixedKern();
  else if(type=="bias")
    pkern = new CBiasKern();

  else if(type=="rbf")
    pkern = new CRbfKern();
  else if(type=="ratquad")
    pkern = new CRatQuadKern();
  else if(type=="matern32")
    pkern = new CMatern32Kern();
  else if(type=="matern52")
    pkern = new CMatern52Kern();

  //else if(type="rbfperiodic")
  //  pkern = new CRbfPeriodicKern();
  //else if(type="gibbsperiodic")
  //  pkern = new CGibbsPeriodicKern();
  
  else if(type=="lin")
    pkern = new CLinKern();
  else if(type=="poly")
    pkern = new CPolyKern();
  else if(type=="mlp")
    pkern = new CMlpKern();
  //else if(type="gibbs")
  //  pkern = new CGibbsKern():

  else if(type=="rbfard")
    pkern = new CRbfardKern();
  //else if(type=="ratquadard")
  //  pkern = new CRatQuadardKern();
  //else if(type=="matern32ard")
  //  pkern = new CMatern32ardKern();
  //else if(type=="matern52ard")
  //  pkern = new CMatern52ardKern();

  else if(type=="linard")
    pkern = new CLinardKern();
  else if(type=="polyard")
    pkern = new CPolyardKern();
  else if(type=="mlpard")
    pkern = new CMlpardKern();
  //else if(type="gibbsard")
  //  pkern = new CGibbsardKern():

  else if(type=="cmpnd")
    pkern = new CCmpndKern();
  else if(type=="tensor")
    pkern = new CTensorKern();
  else
    throw ndlexceptions::StreamFormatError("type", "Unknown kernel type " + type);
  pkern->readParamsFromStream(in);  
  return pkern;
}
void CKern::readParamsFromStream(istream& in)
{
//   string tbaseType = getBaseTypeStream(in);
//   if(tbaseType != getBaseType())
//     throw ndlexceptions::StreamFormatError("baseType", "Error mismatch between saved base type, " + tbaseType + ", and Class base type, " + getType() + ".");
//   string ttype = getTypeStream(in);
//   if(ttype != getType())
//     throw ndlexceptions::StreamFormatError("type", "Error mismatch between saved type, " + ttype + ", and Class type, " + getType() + ".");
  setInputDim(readIntFromStream(in, "inputDim"));
  unsigned int nPars = readIntFromStream(in, "numParams");
  CMatrix par(1, nPars);
  par.fromStream(in);
  if(nPars==getNumParams())
    setParams(par);
  else
    throw ndlexceptions::StreamFormatError("numParams", "Listed number of parameters does not match computed number of parameters.");
  unsigned int numPriors = readIntFromStream(in, "numPriors");
  readPriorsFromStream(in, numPriors);
}
#ifdef _NDLMATLAB
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
  if(mxType!=getType())
  {
    throw ndlexceptions::FileReadError("Error mismatch between saved type, " + mxType + ", and Class type, " + type + ".");
  }
  mxArray* transformArray = mxArrayExtractMxArrayField(matlabArray, "transforms");
  // transforms field.
  transformsFromMxArray(transformArray);
  extractParamFromMxArray(matlabArray);
}
void CKern::extractParamFromMxArray(const mxArray* matlabArray)
{
  inputDim = mxArrayExtractIntField(matlabArray, "inputDimension");
  nParams = mxArrayExtractIntField(matlabArray, "nParams");
  string pName;
  for(unsigned int i=0; i<nParams; i++)
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
  for(unsigned int i=0; i<nParams; i++)
  {      
    pName = getParamName(i);
    mxAddField(matlabArray, pName.c_str());      
    mxSetField(matlabArray, 0, pName.c_str(), convertMxArray(getParam(i))); 
  }
}
void CArdKern::extractParamFromMxArray(const mxArray* matlabArray)
{
  setInputDim(mxArrayExtractIntField(matlabArray, "inputDimension"));
  DIMENSIONMATCH(nParams == mxArrayExtractIntField(matlabArray, "nParams"));
  string pName;
  for(unsigned int i=0; i<nParams; i++)
  {
    pName=getParamName(i);
    if(pName!="inputScale")
      setParam(mxArrayExtractDoubleField(matlabArray, pName), i);  
  }
  scales.fromMxArray(mxArrayExtractMxArrayField(matlabArray, "inputScales"));
}
void CArdKern::addParamToMxArray(mxArray* matlabArray) const
{
  mxAddField(matlabArray, "nParams");
  mxSetField(matlabArray, 0, "nParams", convertMxArray((double)nParams));
  string pName;
  for(unsigned int i=0; i<nParams; i++)
  {      
    pName = getParamName(i);
    if(pName!="inputScale")
    {
      mxAddField(matlabArray, pName.c_str());      
      mxSetField(matlabArray, 0, pName.c_str(), convertMxArray(getParam(i))); 
    }
  }
  mxAddField(matlabArray, "inputScales");
  mxSetField(matlabArray, 0, "inputScales", scales.toMxArray());
}
void CComponentKern::extractParamFromMxArray(const mxArray* matlabArray)
{
  nParams=0;
  mxArray* compArray=mxArrayExtractMxArrayField(matlabArray, "comp");
  unsigned int numComps = mxGetNumberOfElements(compArray);
  mxArray* compElement;
  string kernType;
  unsigned int kernNum = 0;
  for(unsigned int i=0; i<numComps; i++)
  {
    // TODO this is a pretty ugly way of doing things ... perhaps dynamic casting can be used?
    compElement = mxGetCell(compArray, i);
    kernType = mxArrayExtractStringField(compElement, "type");
    if(kernType == "lin")
      kernNum = addKern(new CLinKern(getInputDim()));
    else if(kernType == "rbf")  
      kernNum = addKern(new CRbfKern(getInputDim()));
    else if(kernType == "ratquad")  
      kernNum = addKern(new CRatQuadKern(getInputDim()));
    else if(kernType == "matern32")  
      kernNum = addKern(new CMatern32Kern(getInputDim()));
    else if(kernType == "matern52")  
      kernNum = addKern(new CMatern52Kern(getInputDim()));
    else if(kernType == "mlp")  
      kernNum = addKern(new CMlpKern(getInputDim()));
    else if(kernType == "poly")  
      kernNum = addKern(new CPolyKern(getInputDim()));
    else if(kernType == "polyard")  
      kernNum = addKern(new CPolyardKern(getInputDim()));
    else if(kernType == "linard")  
      kernNum = addKern(new CLinardKern(getInputDim()));
    else if(kernType == "rbfard")  
      kernNum = addKern(new CRbfardKern(getInputDim()));
    else if(kernType == "mlpard")  
      kernNum = addKern(new CMlpardKern(getInputDim()));
    else if(kernType == "white")
      kernNum = addKern(new CWhiteKern(getInputDim()));
    else if(kernType == "bias")
      kernNum = addKern(new CBiasKern(getInputDim()));
    else if(kernType == "cmpnd")
      kernNum = addKern(new CCmpndKern(getInputDim()));
    else if(kernType == "tensor")
      kernNum = addKern(new CTensorKern(getInputDim()));
    else
      throw ndlexceptions::Error("Unknown kernel type.");
    components[kernNum]->fromMxArray(compElement);
  }
}
void CComponentKern::addParamToMxArray(mxArray* matlabArray) const
{
  // Add comp field to mxArray.
  int dims[1];
  mxAddField(matlabArray, "comp");
  dims[0] = components.size();
  mxArray* compArray = mxCreateCellArray(1, dims);
  mxSetField(matlabArray, 0, "comp", 
	     compArray); 
  mxArray* currentKern;
  for(size_t i=0; i<components.size(); i++)
  {
    currentKern = components[i]->toMxArray();
    mxAddField(currentKern, "index");
    mxSetCell(compArray, i, currentKern);
  }
  
}
#endif /* _NDLMATLAB */
