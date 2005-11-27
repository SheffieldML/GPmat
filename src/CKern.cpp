#include "CKern.h"

using namespace std;

/*void CKern::initialiseKern(int inDim) 
{
  setInputDim(inDim);
}
void CKern::initialiseKern(const CMatrix& X)
{
  setInputDim(X.getCols());
}
// Copy constructor
void CKern::initialiseKern(const CKern& kern)
{
  // TODO I think this is messed up ...copy vectors and priors
  // transforms=kern.transforms;
  // TODO  priors=kern.priors;
  inputDim=kern.inputDim;
}*/
ostream& CKern::display(ostream& os) const
{
  os << getName() << " kernel:" << endl;
  for(int i=0; i<nParams; i++)
    {
      os << getParamName(i) << ": " << getParam(i) << endl;
    }
  return os;
}
void CKern::writeParamsToStream(ostream& out) const
{
  out << "numParams=" << getNumParams() << endl;
  out << "numFeatures=" << getInputDim() << endl;
  for(int i=0; i<getNumParams()-1; i++)
    out << getParam(i) << " ";
  out << getParam(getNumParams()-1) << endl;
  out << "numPriors=" << getNumPriors() << endl;
  writePriorsToStream(out);
}
void CKern::getGradPrior(CMatrix& G) const
{
}
void CKern::getPriorLogProb(CMatrix& G) const
{
}

void CKern::getGradTransParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd, bool regularise) const
{
  assert(g.getRows()==1);
  assert(g.getCols()==getNumParams());
  getGradParams(g, X, cvGrd, regularise);
  double val;
  double param;
  for(int i=0; i<getNumTransforms(); i++)
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

// The Compound kernel.
CCmpndKern::CCmpndKern()
{
}
CCmpndKern::CCmpndKern(int inDim)
{
  setInputDim(inDim);
}
CCmpndKern::CCmpndKern(const CMatrix& X) 
{
  setInputDim(X.getCols());
}  
// Class destructor
CCmpndKern::~CCmpndKern()
{
}
CCmpndKern::CCmpndKern(const CCmpndKern& kern) : components(kern.components)
{
  setInputDim(kern.getInputDim());
  for(size_t i=0; i<components.size(); i++)
    addKern(components[i]->clone()); 
  //  components(kern.components);// = kern.components;
}
void CCmpndKern::setInitParam()
{
  setType("cmpnd");
  setName("compound");
  nParams=0;
}
double CCmpndKern::diagComputeElement(const CMatrix& X, int index) const
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
  for(size_t i=0; i < components.size(); i++)
    {
      components[i]->diagCompute(dStore, X);
      d.axpy(dStore, 1.0);
    }
}
void CCmpndKern::setParam(double val, int paramNo)
{
  int start = 0;
  int end = 0;
  for(size_t i=0; i<components.size(); i++)
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
  for(size_t i=0; i<components.size(); i++)
    {
      end = start+components[i]->getNumParams()-1;
      if(paramNo <= end)
	return components[i]->getParam(paramNo-start);
      start = end + 1;
    }
    return -1;
}
string CCmpndKern::getParamName(int paramNo) const
{
  int start = 0;
  int end = 0;
  for(size_t i=0; i<components.size(); i++)
    {
      end = start+components[i]->getNumParams()-1;
      if(paramNo <= end)
	return components[i]->getType() + components[i]->getParamName(paramNo-start);
      start = end + 1;
    }
    return "";
}
void CCmpndKern::getGradX(vector<CMatrix*>& gX, const CMatrix& X, const CMatrix& X2, const bool addG) const
{
  assert(gX.size()==X.getRows());
  assert(X.getCols()==X2.getCols());
  if(!addG)
    for(int i=0; i<gX.size(); i++)
      gX[i]->zeros();
  for(int i=0; i<components.size(); i++)
    components[i]->getGradX(gX, X, X2, true);
}
void CCmpndKern::getDiagGradX(CMatrix& gX, const CMatrix& X, const bool addG) const
{
  assert(gX.dimensionsMatch(X));
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

double CCmpndKern::computeElement(const CMatrix& X1, int index1, 
			   const CMatrix& X2, int index2) const
{
  double y=0.0;
  for(size_t i=0; i<components.size(); i++)
    y+=components[i]->computeElement(X1, index1, X2, index2);
  return y;
}

void CCmpndKern::compute(CMatrix& K, const CMatrix& X) const
{
  assert(K.rowsMatch(X));
  assert(K.isSquare());
  K.zeros();
  CMatrix K2(K.getRows(), K.getCols());
  for(size_t i=0; i<components.size(); i++)
    {
      components[i]->compute(K2, X);
      K.axpy(K2, 1.0);
    }
}

void CCmpndKern::compute(CMatrix& K, const CMatrix& X, const CMatrix& X2) const
{
  assert(K.rowsMatch(X));
  assert(K.getCols()==X2.getRows());
  //CMatrix K2(K.getRows(), K.getCols());
  for(int i=0; i<X.getRows(); i++)
    {
      for(int j=0; j<X2.getRows(); j++)
	{
	  double kval=0.0;
	  for(int k=0; k<components.size(); k++)
	    kval+=components[k]->computeElement(X, i, X2, j);
	  K.setVal(kval, i, j);
	}
    }
}
void CCmpndKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& covGrad, bool regularise) const
{
  assert(g.getRows()==1);
  assert(g.getCols()==getNumParams());
  int start = 0;
  int end = 0;
  for(int i=0; i<components.size(); i++)
    {
      end = start+components[i]->getNumParams()-1;
      CMatrix subg(1, components[i]->getNumParams());
      components[i]->getGradParams(subg, X, covGrad, regularise);
      g.setMatrix(0, start, subg);
      start = end+1;
    }
}
double CCmpndKern::getGradParam(int index, const CMatrix& X, const CMatrix& covGrad) const
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
    return -1;
}
double CCmpndKern::priorLogProb() const
{
  double L = 0.0;
  for(int i=0; i<components.size(); i++)
    {
      L+=components[i]->priorLogProb();
    }
  return L;
}
void CCmpndKern::addPrior(CDist* prior, int index) 
{
  cerr << "Error cannot add priors to compound kernels directly, please add to the components." << endl;
  exit(1);
}

void CCmpndKern::updateX(const CMatrix& X)
{
  for(size_t i=0; i<components.size(); i++)
    components[i]->updateX(X);
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
void CCmpndKern::readParamsFromStream(istream& in) 
{
  string line;
  vector<string> tokens;
  // first line is number of kernels.
  getline(in, line);
  ndlstrutil::tokenise(tokens, line, "=");
  if(tokens.size()!=2 || tokens[0]!="numKerns")
    throw ndlexceptions::FileFormatError();
  int numKerns=atol(tokens[1].c_str());
  tokens.clear();
  getline(in, line);
  ndlstrutil::tokenise(tokens, line, "=");
  if(tokens.size()!=2 || tokens[0]!="numFeatures")
    throw ndlexceptions::FileFormatError();
  int numFeat=atol(tokens[1].c_str());
  setInputDim(numFeat);
  for(int i=0; i<numKerns; i++)
    {
      addKern(readKernFromStream(in));
    }
  
}
void CCmpndKern::writeParamsToStream(ostream& out) const
{
  out << "numKerns=" << components.size() << endl;
  out << "numFeatures=" << getInputDim() << endl;
  for(int i=0; i<components.size(); i++)
    {
      CKern* pkern = components[i];
      writeKernToStream(*pkern, out);
    }
}
// the white noise kernel.
CWhiteKern::CWhiteKern()
{
}
CWhiteKern::CWhiteKern(int inDim)
{
  setInputDim(inDim);
}
CWhiteKern::CWhiteKern(const CMatrix& X)
{
  setInputDim(X.getCols());
}  
// Class destructor
CWhiteKern::~CWhiteKern()
{
}
CWhiteKern::CWhiteKern(const CWhiteKern& kern)
{
  setInputDim(kern.getInputDim());
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
}

inline double CWhiteKern::diagComputeElement(const CMatrix& X, int index) const
{
  return variance;
}
void CWhiteKern::diagCompute(CMatrix& d, const CMatrix& X) const
{
  assert(d.getCols()==1);
  assert(X.rowsMatch(d));
  d.setVals(variance);
}
void CWhiteKern::setParam(double val, int paramNo)
{
  assert(paramNo==0);
  switch(paramNo)
    {
    case 0:
      variance = val;
      break;
    default:
      {
	cerr<<"Parameter doesn't exist.";
	exit(1);
      }
    }
}
// Parameters are kernel parameters
double CWhiteKern::getParam(int paramNo) const
{
  assert(paramNo==0);
  switch(paramNo)
    {
    case 0:
      return variance;
      break;
    default:
      {
	cerr<<"Parameter doesn't exist.";
	exit(1);
      }
    }
}
void CWhiteKern::getGradX(vector<CMatrix*>& gX, const CMatrix& X, const CMatrix& X2, const bool addG) const
{
  assert(gX.size()==X.getRows());
  assert(X.getCols()==X2.getCols());
  if(!addG)
  {
    for(int i=0; i<gX.size(); i++)
	{
	  assert(gX[i]->getRows()==X2.getRows());
	  assert(gX[i]->getCols()==X2.getCols());
	  gX[i]->zeros();
    }
  }
}
void CWhiteKern::getDiagGradX(CMatrix& gX, const CMatrix& X, const bool addG) const
{
  assert(gX.dimensionsMatch(X));
  if(!addG)
    gX.zeros();
}
double CWhiteKern::getWhite() const
{
  return variance;
}

inline double CWhiteKern::computeElement(const CMatrix& X1, int index1,
				  const CMatrix& X2, int index2) const
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
double CWhiteKern::getGradParam(int index, const CMatrix& X, const CMatrix& covGrad) const
{
  assert(index==0);
  return trace(covGrad);
}

// the bias kernel.
CBiasKern::CBiasKern()
{
}
CBiasKern::CBiasKern(int inDim)
{
  setInputDim(inDim);
}
CBiasKern::CBiasKern(const CMatrix& X)
{
  setInputDim(X.getCols());
}  
// Class destructor
CBiasKern::~CBiasKern()
{
}
CBiasKern::CBiasKern(const CBiasKern& kern)
{
  setInputDim(kern.getInputDim());
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
}

double CBiasKern::diagComputeElement(const CMatrix& X, int index) const
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
void CBiasKern::setParam(double val, int paramNo)
{
  assert(paramNo==0);
  switch(paramNo)
    {
    case 0:
      variance = val;
      break;
    default:
      {
	cerr<<"Parameter doesn't exist.";
	exit(1);
      }
    }
}
double CBiasKern::getParam(int paramNo) const
{
  assert(paramNo==0);
  switch(paramNo)
    {
    case 0:
      return variance;
      break;
    default:
      {
	cerr<<"Parameter doesn't exist.";
	exit(1);
      }
    }
}
void CBiasKern::getGradX(vector<CMatrix*>& gX, const CMatrix& X, const CMatrix& X2, const bool addG) const
{
  assert(gX.size()==X.getRows());
  assert(X.getCols()==X2.getCols());
  if(!addG)
  {
    for(int i=0; i<gX.size(); i++)
	{
	  assert(gX[i]->getRows()==X2.getRows());
	  assert(gX[i]->getCols()==X2.getCols());
	  gX[i]->zeros();
	}
  }
}
void CBiasKern::getDiagGradX(CMatrix& gX, const CMatrix& X, const bool addG) const
{
  assert(gX.dimensionsMatch(X));
  if(!addG)
    gX.zeros();
}
double CBiasKern::getWhite() const
{
  return 0.0;
}

inline double CBiasKern::computeElement(const CMatrix& X1, int index1, 
				 const CMatrix& X2, int index2) const
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
double CBiasKern::getGradParam(int index, const CMatrix& X, const CMatrix& covGrad) const 
{
  assert(index==0);
  return sum(covGrad);
}
// the RBF kernel.
CRbfKern::CRbfKern() : updateXused(false)
{
}
CRbfKern::CRbfKern(int inDim) : updateXused(false)
{
  setInputDim(inDim);
}
CRbfKern::CRbfKern(const CMatrix& X) : updateXused(false)
{
  setInputDim(X.getCols());
}  
// Class destructor
CRbfKern::~CRbfKern()
{
}
CRbfKern::CRbfKern(const CRbfKern& kern)
{
  setInputDim(kern.getInputDim());
  variance = kern.variance;
  inverseWidth = kern.inverseWidth;
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
}

inline double CRbfKern::diagComputeElement(const CMatrix& X, int index) const
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
void CRbfKern::setParam(double val, int paramNo)
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
      {
	cerr << "Parameter doesn't exist.";
	exit(1);
      }
    }
}
double CRbfKern::getParam(int paramNo) const
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
      {
	cerr << "Parameter doesn't exist.";
	exit(1);
      }
    }
}

void CRbfKern::getGradX(vector<CMatrix*>& gX, const CMatrix& X, const CMatrix& X2, const bool addG) const
{
  // WVB: always called with X=X2, so I'm not sure why both args are needed.
  // I guess the idea is that the two x vals _could_ be rows living in different
  // physical matrices
  
  // NDL: Mostly called with X=X2, but can also be called with X!=X2
  // when a submatrix of the kernel matrix is involved.

  // This computes a rank 3 tensor with 
  //   gX(i,k,j) = d kern(x_i,x_k)/ d x_component_j
  //
  // The full gradient of dK/dX is a rank 4 tensor, but the row in the 
  // denominator justs adds a [delta(i,j) or delta(k,j)] factor:
  //   gX(i,k,j,l) = gX(i,k,j)*[delta(i,j) or delta(k,j)]
  // I.e. the value is zero unless i or k is equal to j.

  // There is an underlying assumption in getGradX that 
  // kernel entry i,j only depends on row i of X and row j of X2.

  assert(gX.size()==X.getRows());
  assert(X.getCols()==X2.getCols());
  double wi2= 0.5*inverseWidth;
  double pf = variance*inverseWidth;
  for(int i=0; i<X.getRows(); i++)
  {
    assert(gX[i]->getRows()==X2.getRows());
    assert(gX[i]->getCols()==X2.getCols());
    for(int k=0; k<X2.getRows(); k++)
	{
	  double n2 = X.dist2Row(i, X2, k);
	  for(int j=0; j<X2.getCols(); j++)
      {
        double val = pf*(X2.getVal(k, j)-X.getVal(i, j))*exp(-n2*wi2);
        if(addG)
          gX[i]->addVal(val, k, j);
        else 
          gX[i]->setVal(val, k, j);
      }
	}
  }
}
void CRbfKern::getDiagGradX(CMatrix& gX, const CMatrix& X, const bool addG) const
{
  assert(gX.dimensionsMatch(X));
  if(!addG)
    gX.zeros();
}
double CRbfKern::getWhite() const
{
  return 0.0;
}

double CRbfKern::computeElement(const CMatrix& X1, int index1, 
			 const CMatrix& X2, int index2) const
{
  double k = X1.dist2Row(index1, X2, index2);
  k = 0.5*k*inverseWidth;
  k = variance*exp(-k);
  return k;
}

void CRbfKern::updateX(const CMatrix& X)
{
  updateXused = true;
  Xdists.resize(X.getRows(),X.getRows());
  double halfInverseWidth=0.5*inverseWidth;
  int nrows = X.getRows();
  for(int j=0; j<nrows; j++)
  {
    Xdists.setVal(0,j,j);
    for(int i=0; i<j; i++)
    {
      double dist2 = X.dist2Row(i, X, j);
      Xdists.setVal(dist2,i,j);
      Xdists.setVal(exp(-dist2*halfInverseWidth),j,i);
    }
  }
}

void CRbfKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& covGrad, bool regularise) const
{
  assert(g.getRows()==1);
  assert(g.getCols()==nParams);
  assert(covGrad.isSymmetric());
  double g1=0.0;
  double g2=0.0;
  double halfInverseWidth=0.5*inverseWidth;
  double halfVariance=0.5*variance;

  int nrows = X.getRows();
  for(int j=0; j<nrows; j++)
  {
    g2 += covGrad.getVal(j,j);
    for(int i=0; i<j; i++)
    {
      double k = 0;
      double dist2 = 0;
      if(updateXused) // WVB's mod for precomputing parts of the kernel.
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

double CRbfKern::getGradParam(int index, const CMatrix& X, const CMatrix& covGrad) const
{
  assert(index>=0);
  assert(index<nParams);
  cerr << "Error getGradParam is not currently implemented for CRbfKern" << endl;
  exit(1);
}


// the Linear kernel.
CLinKern::CLinKern()
{
}
CLinKern::CLinKern(int inDim)
{
  setInputDim(inDim);
}
CLinKern::CLinKern(const CMatrix& X)
{
  setInputDim(X.getCols());
}  
// Class destructor
CLinKern::~CLinKern()
{
}
CLinKern::CLinKern(const CLinKern& kern)
{
  setInputDim(kern.getInputDim());
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
}

double CLinKern::diagComputeElement(const CMatrix& X, int index1) const
{
  return variance*X.norm2Row(index1);  
}
// Parameters are kernel parameters
void CLinKern::setParam(double val, int paramNo)
{
  assert(paramNo==0);
  switch(paramNo)
    {
    case 0:
      variance = val;
      break;
    default:
      {
	cerr << "Parameter doesn't exist.";
	exit(1);
      }
    }
}
double CLinKern::getParam(int paramNo) const
{
  assert(paramNo==0);
  switch(paramNo)
    {
    case 0:
      return variance;
      break;
    default:
      {
	cerr << "Parameter doesn't exist.";
	exit(1);
      }
    }
}
void CLinKern::getGradX(vector<CMatrix*>& gX, const CMatrix& X, const CMatrix& X2, const bool addG) const
{
  assert(gX.size()==X.getRows());
  assert(X.getCols()==X2.getCols());
  for(int i=0; i<X.getRows(); i++)
    {
      assert(gX[i]->getRows()==X2.getRows());
      assert(gX[i]->getCols()==X2.getCols());
      for(int k=0; k<X2.getRows(); k++)
	{
	  for(int j=0; j<X2.getCols(); j++)
	    {
	      double val = 0.0;
	      if(addG)
		val = gX[i]->getVal(k, j);
	      val += variance*X2.getVal(k, j);
	      gX[i]->setVal(val, k, j);
	    }
	}
    }
}
void CLinKern::getDiagGradX(CMatrix& gX, const CMatrix& X, const bool addG) const
{
  assert(gX.dimensionsMatch(X));
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

double CLinKern::computeElement(const CMatrix& X1, int index1, 
			  const CMatrix& X2, int index2) const
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
double CLinKern::getGradParam(int index, const CMatrix& X, const CMatrix& covGrad) const
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
  

CMlpKern::CMlpKern()
{
}
CMlpKern::CMlpKern(int inDim)
{
  setInputDim(inDim);
}
CMlpKern::CMlpKern(const CMatrix& X)
{
  setInputDim(X.getCols());
}  
// Class destructor
CMlpKern::~CMlpKern()
{
}
CMlpKern::CMlpKern(const CMlpKern& kern)
{
  setInputDim(kern.getInputDim());
  variance = kern.variance;
  weightVariance = kern.weightVariance;
  biasVariance = kern.biasVariance;
}
void CMlpKern::setInitParam()
{
  nParams = 3;
  setType("mlp");
  setName("MLP");
  setParamName("weightVariance", 0);
  weightVariance = 10.0;
  setParamName("biasVariance", 1);
  biasVariance = 10.0;
  setParamName("variance", 2);
  variance = 1.0;
  addTransform(new CNegLogLogitTransform, 0);
  addTransform(new CNegLogLogitTransform, 1);
  addTransform(new CNegLogLogitTransform, 2);
}

inline double CMlpKern::diagComputeElement(const CMatrix& X, int index) const
{
  double numer=weightVariance*X.norm2Row(index)+biasVariance;  
  double denom = numer+1.0;
  return variance*asin(numer/denom);
}
// Parameters are kernel parameters
void CMlpKern::setParam(double val, int paramNo)
{
  assert(paramNo < nParams);
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
      {
	cerr << "Parameter doesn't exist.";
	exit(1);
      }
    }
}
double CMlpKern::getParam(int paramNo) const
{
  assert(paramNo < nParams);
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
      {
	cerr << "Parameter doesn't exist.";
	exit(1);
      }
    }
}

void CMlpKern::getGradX(vector<CMatrix*>& gX, const CMatrix& X, const CMatrix& X2, const bool addG) const
{
  assert(gX.size()==X.getRows());
  assert(X.getCols()==X2.getCols());
  for(int i=0; i<X.getRows(); i++)
    {
      assert(gX[i]->getRows()==X2.getRows());
      assert(gX[i]->getCols()==X2.getCols());
      double denomPart1 = weightVariance*X.norm2Row(i) + biasVariance + 1.0;
      for(int k=0; k<X2.getRows(); k++)
	{
	  double numer=weightVariance*X.dotRowRow(i, X2, k) + biasVariance;
	  double denomPart2=weightVariance*X2.norm2Row(k) + biasVariance + 1.0;
	  double denom=sqrt(denomPart1*denomPart2);
	  double arg = numer/denom;
	  double kval=variance*weightVariance/sqrt(1-arg*arg);
	  
	  for(int j=0; j<X2.getCols(); j++)
	    {
	      double val = 0.0;
	      if(addG)
		val = gX[i]->getVal(k, j);
	      val += (X2.getVal(k, j)/denom-X.getVal(i, j)*denomPart2*numer/(denom*denom*denom))*kval;
	      gX[i]->setVal(val, k, j);
	    }
	}
    }

}
void CMlpKern::getDiagGradX(CMatrix& gX, const CMatrix& X, const bool addG) const
{
  assert(gX.dimensionsMatch(X));
  for(int i=0; i<X.getRows(); i++)
    {
      double numer=weightVariance*X.norm2Row(i) + biasVariance;
      double denom=numer+1.0;
      double arg=numer/denom;
      double kval=variance*weightVariance/sqrt(1-arg*arg);
      for(int j=0; j<X.getCols(); j++)
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

double CMlpKern::computeElement(const CMatrix& X1, int index1, 
			 const CMatrix& X2, int index2) const
{
  double numer= weightVariance*X1.dotRowRow(index1, X2, index2) + biasVariance;
  double denom1=weightVariance*X1.norm2Row(index1)+biasVariance+1.0;  
  double denom2=weightVariance*X2.norm2Row(index2)+biasVariance+1.0;  
  return variance*asin(numer/sqrt(denom1*denom2));
}
void CMlpKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& covGrad, bool regularise) const
{
  assert(g.getRows()==1);
  assert(g.getCols()==nParams);
  double g1=0.0;
  double g2=0.0;
  double g3=0.0;
  innerProd.resize(1, X.getRows());

  // do off diagonal gradients first.
  for(int i=0; i<X.getRows(); i++)
    {
      innerProd.setVal(X.norm2Row(i), i);
      for(int j=0; j<i; j++)
	{	  
	  double crossProd=X.dotRowRow(i, X, j);
	  double numer=weightVariance*crossProd + biasVariance;
	  double denomi=weightVariance*innerProd.getVal(i)+biasVariance+1.0;  
	  double denomj=weightVariance*innerProd.getVal(j)+biasVariance+1.0;  
	  double denom = sqrt(denomi*denomj);
	  double arg = numer/denom;
	  double baseCovGrad = variance/sqrt(1-arg*arg)*covGrad.getVal(i, j);
	  double denom3=denom*denom*denom;
	  g1+=baseCovGrad*(crossProd/denom-0.5*numer/denom3*(denomi*innerProd.getVal(j) + innerProd.getVal(i)*denomj));
	  g2+=baseCovGrad*(1.0/denom-0.5*numer/denom3*((innerProd.getVal(i)+innerProd.getVal(j))*weightVariance +2.0*biasVariance+2.0));
	  g3+=asin(arg)*covGrad.getVal(i, j);
	}
    }
  // double the result due to symmetry.
  g1*=2;
  g2*=2;
  g3*=2;
  // add effect of diagonals.
  for(int i=0; i<X.getRows(); i++)
    {
      double numer = weightVariance*innerProd.getVal(i)+biasVariance;
      double denom = numer+1.0;
      double denom3=denom*denom*denom;
      double arg = numer/denom;
      double baseCovGrad = variance/sqrt(1-arg*arg)*covGrad.getVal(i, i);
      g1+=baseCovGrad*(innerProd.getVal(i)/denom-0.5*numer/denom3*(2*denom*innerProd.getVal(i)));
      
      g2+=baseCovGrad*(1.0/denom-0.5*numer/denom3*(2.0*weightVariance*innerProd.getVal(i)+2.0*biasVariance+2.0));
      g3+=asin(arg)*covGrad.getVal(i, i);

    }
  g.setVal(g1, 0);
  g.setVal(g2, 1);
  g.setVal(g3, 2);
  if(regularise)
    addPriorGrad(g);

}
double CMlpKern::getGradParam(int index, const CMatrix& X, const CMatrix& covGrad) const
{
  assert(index>=0);
  assert(index<nParams);
  cerr << "Error getGradParam is not currently implemented for CMlpKern" << endl;
  exit(1);
}

CPolyKern::CPolyKern()
{
}
CPolyKern::CPolyKern(int inDim)
{
  setInputDim(inDim);
}
CPolyKern::CPolyKern(const CMatrix& X)
{
  setInputDim(X.getCols());
}  
// Class destructor
CPolyKern::~CPolyKern()
{
}
CPolyKern::CPolyKern(const CPolyKern& kern)
{
  setInputDim(kern.getInputDim());
  variance = kern.variance;
  weightVariance = kern.weightVariance;
  biasVariance = kern.biasVariance;
  degree = kern.degree;
}
void CPolyKern::writeParamsToStream(ostream& out) const
{
  out << "numParams=" << getNumParams() << endl;
  out << "numFeatures=" << getInputDim() << endl;
  double deg = getDegree();
  if((deg - (int)deg)==0)
    out << "degree=" << (int)deg << endl;
  else
    out << "degree=" << deg << endl;
  for(int i=0; i<getNumParams()-1; i++)
    out << getParam(i) << " ";
  out << getParam(getNumParams()-1) << endl;
}
void CPolyKern::readParamsFromStream(istream& in)
{
  string line;
  vector<string> tokens;
  getline(in, line);
  ndlstrutil::tokenise(tokens, line, "=");
  if(tokens.size()>2 || tokens[0]!="numParams")
    throw ndlexceptions::FileFormatError();
  int numParams=atol(tokens[1].c_str());
  tokens.clear();
  getline(in, line);
  ndlstrutil::tokenise(tokens, line, "=");
  if(tokens.size()>2 || tokens[0]!="numFeatures")
    throw ndlexceptions::FileFormatError();
  int numFeatures=atol(tokens[1].c_str());
  tokens.clear();
  getline(in, line);
  ndlstrutil::tokenise(tokens, line, "=");
  if(tokens.size()>2 || tokens[0]!="degree")
    throw ndlexceptions::FileFormatError();
  setDegree(atol(tokens[1].c_str()));
  CMatrix par(1, numParams);
  tokens.clear();
  getline(in, line);
  ndlstrutil::tokenise(tokens, line, " ");
  for(int i=0; i<numParams; i++)
    par.setVal(atof(tokens[i].c_str()), i);
  setInputDim(numFeatures);
  if(numParams==getNumParams())
    setParams(par);
  else
    throw ndlexceptions::FileFormatError();
}
void CPolyKern::setInitParam()
{
  nParams = 3;
  setType("poly");
  setName("Polynomial");
  setParamName("weightVariance", 0);
  weightVariance = 1.0;
  setParamName("biasVariance", 1);
  biasVariance = 1.0;
  setParamName("variance", 2);
  variance = 1.0;
  addTransform(new CNegLogLogitTransform, 0);
  addTransform(new CNegLogLogitTransform, 1);
  addTransform(new CNegLogLogitTransform, 2);
  degree = 2.0;
}

inline double CPolyKern::diagComputeElement(const CMatrix& X, int index) const
{
  double arg=weightVariance*X.norm2Row(index)+biasVariance;  
  return variance*pow(arg, degree);
}
// Parameters are kernel parameters
void CPolyKern::setParam(double val, int paramNo)
{
  assert(paramNo < nParams);
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
      {
	cerr << "Parameter doesn't exist.";
	exit(1);
      }
    }
}
double CPolyKern::getParam(int paramNo) const
{
  assert(paramNo < nParams);
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
      {
	cerr << "Parameter doesn't exist.";
	exit(1);
      }
    }
}

void CPolyKern::getGradX(vector<CMatrix*>& gX, const CMatrix& X, const CMatrix& X2, const bool addG) const
{
  assert(gX.size()==X.getRows());
  assert(X.getCols()==X2.getCols());
  for(int i=0; i<X.getRows(); i++)
    {
      assert(gX[i]->getRows()==X2.getRows());
      assert(gX[i]->getCols()==X2.getCols());
      for(int k=0; k<X2.getRows(); k++)
	{
	  double arg=weightVariance*X.dotRowRow(i, X2, k) + biasVariance;
	  double kval=degree*variance*weightVariance*pow(arg, degree-1);
	  for(int j=0; j<X2.getCols(); j++)
	    {
	      double val = 0.0;
	      if(addG)
		val = gX[i]->getVal(k, j);
	      val += kval*X2.getVal(k, j);
	      gX[i]->setVal(val, k, j);
	    }
	}
    }
}
void CPolyKern::getDiagGradX(CMatrix& gX, const CMatrix& X, const bool addG) const
{
  assert(gX.dimensionsMatch(X));
  for(int i=0; i<X.getRows(); i++)
    {
      double arg=weightVariance*X.norm2Row(i) + biasVariance;
      double kval=degree*variance*weightVariance*pow(arg, degree-1);
      for(int j=0; j<X.getCols(); j++)
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

double CPolyKern::computeElement(const CMatrix& X1, int index1, 
			 const CMatrix& X2, int index2) const
{
  double arg=weightVariance*X1.dotRowRow(index1, X2, index2) + biasVariance;
  return variance*pow(arg, degree);
}
void CPolyKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& covGrad, bool regularise) const
{
  assert(g.getRows()==1);
  assert(g.getCols()==nParams);
  double g1=0.0;
  double g2=0.0;
  double g3=0.0;
  innerProd.resize(1, X.getRows());
  // do off diagonal gradients first.
  for(int i=0; i<X.getRows(); i++)
    {
      innerProd.setVal(X.norm2Row(i),i);
      for(int j=0; j<i; j++)
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
  for(int i=0; i<X.getRows(); i++)
    {
      double arg = weightVariance*innerProd.getVal(i)+biasVariance;
      double base = variance*degree*pow(arg, degree-1)*covGrad.getVal(i, i);	  
      g1+=innerProd.getVal(i)*base;
      g2+=base;
      g3+=pow(arg,degree)*covGrad.getVal(i, i);

    }
  g.setVal(g1, 0);
  g.setVal(g2, 1);
  g.setVal(g3, 2);
  if(regularise)
    addPriorGrad(g);

}
double CPolyKern::getGradParam(int index, const CMatrix& X, const CMatrix& covGrad) const
{
  assert(index>=0);
  assert(index<nParams);
  cerr << "Error getGradParam is not currently implemented for CPolyKern" << endl;
  exit(1);
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
CLinardKern::CLinardKern()
{
}
CLinardKern::CLinardKern(int inDim)
{
  setInputDim(inDim);
}
CLinardKern::CLinardKern(const CMatrix& X)
{
  setInputDim(X.getCols());
}  
// Class destructor
CLinardKern::~CLinardKern()
{
}
CLinardKern::CLinardKern(const CLinardKern& kern)
{
  setInputDim(kern.getInputDim());
  variance = kern.variance;
  scales = kern.scales;
}
void CLinardKern::setInitParam()
{
  nParams = 1+getInputDim();
  setType("linard");
  setName("linear ARD");
  setParamName("variance", 0);
  variance = 1.0;
  addTransform(new CNegLogLogitTransform, 0);
  scales.resize(1, getInputDim());
  scales.setVals(0.5);
  for(int i=1; i<getInputDim()+1; i++)
    {
      addTransform(new CSigmoidTransform, i);
      string name = "inputScale";
      setParamName(name, i);
    }
}

double CLinardKern::diagComputeElement(const CMatrix& X, int index1) const
{
  double val = 0.0;
  for(int i=0; i<getInputDim(); i++)
    {
      double x=X.getVal(index1, i);
      val += x*x*scales.getVal(i);
    }
  return val*variance;
}
// Parameters are kernel parameters
void CLinardKern::setParam(double val, int paramNo)
{
  assert(paramNo>=0);
  assert(paramNo<nParams);
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
	  cerr << "Parameter doesn't exist.";
	  exit(1);
	}
      
    }
}
double CLinardKern::getParam(int paramNo) const
{
  assert(paramNo>=0);
  assert(paramNo<nParams);
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
	  cerr << "Parameter doesn't exist.";
	  exit(1);
	}
    }
}
void CLinardKern::getGradX(vector<CMatrix*>& gX, const CMatrix& X, const CMatrix& X2, const bool addG) const
{
  assert(gX.size()==X.getRows());
  assert(X.getCols()==X2.getCols());
  for(int i=0; i<X.getRows(); i++)
    {
      assert(gX[i]->getRows()==X2.getRows());
      assert(gX[i]->getCols()==X2.getCols());
      for(int k=0; k<X2.getRows(); k++)
	{
	  for(int j=0; j<X2.getCols(); j++)
	    {
	      double val = 0.0;
	      if(addG)
		val = gX[i]->getVal(k, j);
	      val += variance*X2.getVal(k, j)*scales.getVal(j);
	      gX[i]->setVal(val, k, j);
	    }
	}
    }
}
void CLinardKern::getDiagGradX(CMatrix& gX, const CMatrix& X, const bool addG) const
{
  assert(gX.dimensionsMatch(X));
  if(!addG)
    {
      gX.deepCopy(X);
      double baseScale = 2.0*variance;
      for(int j=0; j<gX.getCols(); j++)
	gX.scaleCol(j, baseScale*scales.getVal(j));
    }
  else
    {
      double baseScale = 2.0*variance;
      for(int j=0; j<gX.getCols(); j++)
	gX.axpyColCol(j, X, j, baseScale*scales.getVal(j));
    }
}
double CLinardKern::getWhite() const
{
  return 0.0;
}

double CLinardKern::computeElement(const CMatrix& X1, int index1, 
			  const CMatrix& X2, int index2) const
{
  double val = 0.0;
  for(int i=0; i<getInputDim(); i++)
    {
      val += X1.getVal(index1, i)*X2.getVal(index2, i)*scales.getVal(i);
    }
  return val*variance;
}

void CLinardKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& covGrad, bool regularise) const
{
  assert(g.getRows()==1);
  assert(g.getCols()==nParams);
  double g1=0.0;
  for(int i=0; i<X.getRows(); i++)
    {
      for(int j=0; j<i; j++)
	{
	  double val=0.0;
	  for(int k=0; k<getInputDim(); k++)
	    val += X.getVal(i, k)*X.getVal(j, k)*scales.getVal(k);
	  g1+=covGrad.getVal(i, j)*val;
	}
    }
  g1*=2.0;
  for(int i=0; i<X.getRows(); i++)
    {
      double val=0.0;
      for(int k=0; k<getInputDim(); k++)
	val += X.getVal(i, k)*X.getVal(i, k)*scales.getVal(k);
      g1+=covGrad.getVal(i, i)*val;
    }
  g.setVal(g1, 0);
  for(int k=0; k<getInputDim(); k++)
    {
      double g2=0.0;
      for(int i=0; i<X.getRows(); i++)
	{
	  for(int j=0; j<i; j++)
	    {
	      g2+=X.getVal(i, k)*X.getVal(j, k)*covGrad.getVal(i, j);
	    }
	}
      g2*=2;
      for(int i=0; i<X.getRows(); i++)
	g2+=X.getVal(i, k)*X.getVal(i, k)*covGrad.getVal(i, i);
      g2*=variance;
      g.setVal(g2, k+1);
    }
  if(regularise)
    addPriorGrad(g);

}
double CLinardKern::getGradParam(int index, const CMatrix& X, const CMatrix& covGrad) const
{
  assert(index>=0);
  assert(index<nParams);
  cerr << "Error getGradParam is not currently implemented for CLinardKern" << endl;
  exit(1);
}

// the RBF ARD kernel.
CRbfardKern::CRbfardKern()
{
}
CRbfardKern::CRbfardKern(int inDim)
{
  setInputDim(inDim);
}
CRbfardKern::CRbfardKern(const CMatrix& X)
{
  setInputDim(X.getCols());
}  
// Class destructor
CRbfardKern::~CRbfardKern()
{
}
CRbfardKern::CRbfardKern(const CRbfardKern& kern)
{
  setInputDim(kern.getInputDim());
  variance = kern.variance;
  inverseWidth = kern.inverseWidth;
  scales = kern.scales;
}
void CRbfardKern::setInitParam()
{
  nParams = 2+getInputDim();
  setType("rbfard");
  setName("RBF ARD");
  setParamName("inverseWidth", 0);
  inverseWidth=1.0;
  addTransform(new CNegLogLogitTransform, 0);
  setParamName("variance", 1);
  variance = 1.0;
  addTransform(new CNegLogLogitTransform, 1);
  scales.resize(1, getInputDim());
  gscales.resize(1, getInputDim());
  scales.setVals(0.5);
  for(int i=2; i<getInputDim()+2; i++)
    {
      addTransform(new CSigmoidTransform, i);
      string name = "inputScale";
      setParamName(name, i);
    }
}

double CRbfardKern::diagComputeElement(const CMatrix& X, int index1) const
{
  return variance;
}
// Parameters are kernel parameters
void CRbfardKern::setParam(double val, int paramNo)
{
  assert(paramNo>=0);
  assert(paramNo<nParams);
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
	  cerr << "Parameter doesn't exist.";
	  exit(1);
	}
      
    }
}
double CRbfardKern::getParam(int paramNo) const
{
  assert(paramNo>=0);
  assert(paramNo<nParams);
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
	  cerr << "Parameter doesn't exist.";
	  exit(1);
	}
    }
}
void CRbfardKern::getGradX(vector<CMatrix*>& gX, const CMatrix& X, const CMatrix& X2, const bool addG) const
{
  assert(gX.size()==X.getRows());
  assert(X.getCols()==X2.getCols());
  double wi2= 0.5*inverseWidth;
  double pf = variance*inverseWidth;
  for(int i=0; i<X.getRows(); i++)
    {
      assert(gX[i]->getRows()==X2.getRows());
      assert(gX[i]->getCols()==X2.getCols());
      for(int k=0; k<X2.getRows(); k++)
	{
	  double n2=0.0;
	  for(int j=0; j<getInputDim(); j++)
	    {
	      double x = X.getVal(i, j);
	      x = x-X2.getVal(k, j);
	      n2+=x*scales.getVal(j)*x;
	    }
	  for(int j=0; j<X2.getCols(); j++)
	    {
	      double val = 0.0;
	      if(addG)
		val = gX[i]->getVal(k, j);
	      val += pf*(X2.getVal(k, j)-X.getVal(i, j))*exp(-n2*wi2)*scales.getVal(j);
	      gX[i]->setVal(val, k, j);	      
	    }
	}
    }
}
void CRbfardKern::getDiagGradX(CMatrix& gX, const CMatrix& X, const bool addG) const
{
  assert(gX.dimensionsMatch(X));
  if(!addG)
    gX.zeros();
}
double CRbfardKern::getWhite() const
{
  return 0.0;
}

double CRbfardKern::computeElement(const CMatrix& X1, int index1, 
			  const CMatrix& X2, int index2) const
{
  double val = 0.0;
  for(int i=0; i<getInputDim(); i++)
    {
      double x = X1.getVal(index1, i);
      x = x-X2.getVal(index2, i);
      val+=x*scales.getVal(i)*x;
    }
  return variance*exp(-val*inverseWidth*0.5);
}

void CRbfardKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& covGrad, bool regularise) const
{
  assert(g.getRows()==1);
  assert(g.getCols()==nParams);
  double g1=0.0;
  double g2=0.0;
  gscales.zeros();
  double halfInverseWidth = 0.5*inverseWidth;
  for(int i=0; i<X.getRows(); i++)
    {
      for(int j=0; j<i; j++)
	{
	  double val = 0.0;
	  for(int k=0; k<getInputDim(); k++)
	    {
	      double x = X.getVal(i, k);
	      x-=X.getVal(j, k);
	      val+=x*scales.getVal(k)*x;
	    }
	  double kCovGrad = exp(-halfInverseWidth*val)*covGrad.getVal(i, j);
	  g1-=0.5*val*kCovGrad*variance;
	  g2+=kCovGrad;
	  for(int k=0; k<getInputDim(); k++)
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
  for(int i=0; i<X.getRows(); i++)
    g2+=covGrad.getVal(i, i);
  gscales.scale(2.0);
  g.setVal(g1, 0);
  g.setVal(g2, 1);
  for(int k=0; k<getInputDim(); k++)
    g.setVal(gscales.getVal(k), k+2);
  if(regularise)
    addPriorGrad(g);

}
double CRbfardKern::getGradParam(int index, const CMatrix& X, const CMatrix& covGrad) const
{
  assert(index>=0);
  assert(index<nParams);
  cerr << "Error getGradParam is not currently implemented for CRbfardKern" << endl;
  exit(1);
}

// the MLP ARD kernel.
CMlpardKern::CMlpardKern()
{
}
CMlpardKern::CMlpardKern(int inDim)
{
  setInputDim(inDim);
}
CMlpardKern::CMlpardKern(const CMatrix& X)
{
  setInputDim(X.getCols());
}  
// Class destructor
CMlpardKern::~CMlpardKern()
{
}
CMlpardKern::CMlpardKern(const CMlpardKern& kern)
{
  setInputDim(kern.getInputDim());
  variance = kern.variance;
  weightVariance=kern.weightVariance;
  biasVariance=kern.biasVariance;
  scales = kern.scales;
}
void CMlpardKern::setInitParam()
{
  nParams = 3+getInputDim();
  setType("mlpard");
  setName("MLP ARD");
  setParamName("weightVariance", 0);
  weightVariance=10.0;
  addTransform(new CNegLogLogitTransform, 0);
  setParamName("biasVariance", 1);
  biasVariance=10.0;
  addTransform(new CNegLogLogitTransform, 1);
  setParamName("variance", 2);
  variance = 1.0;
  addTransform(new CNegLogLogitTransform, 2);
  scales.resize(1, getInputDim());
  gscales.resize(1, getInputDim());
  scales.setVals(0.5);
  for(int i=3; i<getInputDim()+3; i++)
    {
      addTransform(new CSigmoidTransform, i);
      string name = "inputScale";
      setParamName(name, i);
    }
}

double CMlpardKern::diagComputeElement(const CMatrix& X, int index) const
{
  double val=0.0;
  for(int i=0; i<getInputDim(); i++)
    {
      double x = X.getVal(index, i);
      val+=x*x*scales.getVal(i);
    }
  double numer=weightVariance*val+biasVariance;  
  double denom = numer+1.0;
  return variance*asin(numer/denom);
}
// Parameters are kernel parameters
void CMlpardKern::setParam(double val, int paramNo)
{
  assert(paramNo>=0);
  assert(paramNo<nParams);
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
	  cerr << "Parameter doesn't exist.";
	  exit(1);
	}
      
    }
}
double CMlpardKern::getParam(int paramNo) const
{
  assert(paramNo>=0);
  assert(paramNo<nParams);
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
	  cerr << "Parameter doesn't exist.";
	  exit(1);
	}
    }
}
void CMlpardKern::getGradX(vector<CMatrix*>& gX, const CMatrix& X, const CMatrix& X2, const bool addG) const
{
  assert(gX.size()==X.getRows());
  assert(X.getCols()==X2.getCols());
  for(int i=0; i<X.getRows(); i++)
    {
      assert(gX[i]->getRows()==X2.getRows());
      assert(gX[i]->getCols()==X2.getCols());
      double val=0.0;
      for(int j=0; j<getInputDim(); j++)
	{
	  double x = X.getVal(i, j);
	  val+=x*x*scales.getVal(j);
	}
      double denomPart1=weightVariance*val+biasVariance+1.0;  
      for(int k=0; k<X2.getRows(); k++)
	{
	  double val=0.0;
	  double val2=0.0;
	  for(int j=0; j<getInputDim(); j++)
	    {
	      double xk = X2.getVal(k, j);
	      double xi = X.getVal(i, j);
	      double scalesX = xk*scales.getVal(j);
	      val+=xk*scalesX;
	      val2+=xi*scalesX;
	    }
	  double numer=weightVariance*val2 + biasVariance;
	  double denomPart2=weightVariance*val+biasVariance+1.0;  
	  double denom=sqrt(denomPart1*denomPart2);
	  double arg = numer/denom;
	  double kval=variance*weightVariance/sqrt(1-arg*arg);
	  
	  for(int j=0; j<X2.getCols(); j++)
	    {
	      double val = 0.0;
	      if(addG)
		val = gX[i]->getVal(k, j);
	      val += (X2.getVal(k, j)/denom-X.getVal(i, j)*denomPart2*numer/(denom*denom*denom))*kval*scales.getVal(j);
	      gX[i]->setVal(val, k, j);
	    }
	}
    }

}
void CMlpardKern::getDiagGradX(CMatrix& gX, const CMatrix& X, const bool addG) const
{
  assert(gX.dimensionsMatch(X));
  for(int i=0; i<X.getRows(); i++)
    {
      double val=0.0;
      for(int j=0; j<getInputDim(); j++)
	{
	  double x = X.getVal(i, j);
	  val+=x*x*scales.getVal(j);
	}
      double numer=weightVariance*val+biasVariance;  
      double denom=numer+1.0;
      double arg=numer/denom;
      double kval=variance*weightVariance/sqrt(1-arg*arg);
      for(int j=0; j<X.getCols(); j++)
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

double CMlpardKern::computeElement(const CMatrix& X1, int index1, 
			  const CMatrix& X2, int index2) const
{
  double valij=0.0;
  double valii=0.0;
  double valjj=0.0;
  for(int i=0; i<getInputDim(); i++)
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

void CMlpardKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& covGrad, bool regularise) const
{
  assert(g.getRows()==1);
  assert(g.getCols()==nParams);
  double g1=0.0;
  double g2=0.0;
  double g3=0.0;
  innerProd.resize(1, X.getRows());
  gscales.zeros();
  innerProd.zeros();
  for(int i=0; i<X.getRows(); i++)
    {
      double val=0.0;
      for(int k=0; k<getInputDim(); k++)
	{
	  double x = X.getVal(i, k);
	  val+=x*x*scales.getVal(k);
	}
      innerProd.setVal(val, i);
    }
  
  for(int i=0; i<X.getRows(); i++)
    {
      for(int j=0; j<i; j++)
	{
	  double val=0.0;
	  for(int k=0; k<getInputDim(); k++)
	    {
	      double xi = X.getVal(i, k);
	      double xj = X.getVal(j, k);
	      val+=xi*xj*scales.getVal(k);
	    }
	  double crossProd=val;
	  double numer=weightVariance*val + biasVariance;
	  double denomi=weightVariance*innerProd.getVal(i)+biasVariance+1.0;  
	  double denomj=weightVariance*innerProd.getVal(j)+biasVariance+1.0;  
	  double denom = sqrt(denomi*denomj);
	  double arg = numer/denom;
	  double baseCovGrad = variance/sqrt(1-arg*arg)*covGrad.getVal(i, j);
	  double denom3=denom*denom*denom;
	  g1+=baseCovGrad*(val/denom-0.5*numer/denom3*(denomi*innerProd.getVal(j) + innerProd.getVal(i)*denomj));
	  g2+=baseCovGrad*(1.0/denom-0.5*numer/denom3*((innerProd.getVal(i)+innerProd.getVal(j))*weightVariance+2.0*biasVariance+2.0));
	  g3+=asin(arg)*covGrad.getVal(i, j);


	  for(int k=0; k<getInputDim(); k++)
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
  for(int i=0; i<X.getRows(); i++)
    {
      double numer = weightVariance*innerProd.getVal(i)+biasVariance;
      double denom = numer+1.0;
      double denom3=denom*denom*denom;
      double arg = numer/denom;
      double baseCovGrad = variance/sqrt(1-arg*arg)*covGrad.getVal(i, i);
      g1+=baseCovGrad*(innerProd.getVal(i)/denom-0.5*numer/denom3*(2*denom*innerProd.getVal(i)));
      
      g2+=baseCovGrad*(1.0/denom-0.5*numer/denom3*(2.0*weightVariance*innerProd.getVal(i)+2.0*biasVariance+2.0));
      g3+=asin(arg)*covGrad.getVal(i, i);

      for(int k=0; k<getInputDim(); k++)
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
  for(int k=0; k<getInputDim(); k++)
    g.setVal(gscales.getVal(k), k+3);
  if(regularise)
    addPriorGrad(g);

}
double CMlpardKern::getGradParam(int index, const CMatrix& X, const CMatrix& covGrad) const
{
  assert(index>=0);
  assert(index<nParams);
  cerr << "Error getGradParam is not currently implemented for CMlpardKern" << endl;
  exit(1);
}

// the POLY ARD kernel.
CPolyardKern::CPolyardKern()
{
}
CPolyardKern::CPolyardKern(int inDim)
{
  setInputDim(inDim);
}
CPolyardKern::CPolyardKern(const CMatrix& X)
{
  setInputDim(X.getCols());
}  
// Class destructor
CPolyardKern::~CPolyardKern()
{
}
CPolyardKern::CPolyardKern(const CPolyardKern& kern)
{
  setInputDim(kern.getInputDim());
  variance = kern.variance;
  weightVariance=kern.weightVariance;
  biasVariance=kern.biasVariance;
  scales = kern.scales;
}
void CPolyardKern::writeParamsToStream(ostream& out) const
{
  out << "numParams=" << getNumParams() << endl;
  out << "numFeatures=" << getInputDim() << endl;
  out << "degree=" << getDegree() << endl;
  for(int i=0; i<getNumParams()-1; i++)
    out << getParam(i) << " ";
  out << getParam(getNumParams()-1) << endl;
}
void CPolyardKern::readParamsFromStream(istream& in)
{
  string line;
  vector<string> tokens;
  getline(in, line);
  ndlstrutil::tokenise(tokens, line, "=");
  if(tokens.size()>2 || tokens[0]!="numParams")
    throw ndlexceptions::FileFormatError();
  int numParams=atol(tokens[1].c_str());
  tokens.clear();
  getline(in, line);
  ndlstrutil::tokenise(tokens, line, "=");
  if(tokens.size()>2 || tokens[0]!="numFeatures")
    throw ndlexceptions::FileFormatError();
  int numFeatures=atol(tokens[1].c_str());
  tokens.clear();
  getline(in, line);
  ndlstrutil::tokenise(tokens, line, "=");
  if(tokens.size()>2 || tokens[0]!="degree")
    throw ndlexceptions::FileFormatError();
  int setDegree(atol(tokens[1].c_str()));
  CMatrix par(1, numParams);
  tokens.clear();
  getline(in, line);
  ndlstrutil::tokenise(tokens, line, " ");
  for(int i=0; i<numParams; i++)
    par.setVal(atof(tokens[i].c_str()), i);
  setInputDim(numFeatures);
  if(numParams==getNumParams())
    setParams(par);
  else
    throw ndlexceptions::FileFormatError();
}
void CPolyardKern::setInitParam()
{
  nParams = 3+getInputDim();
  setType("polyard");
  setName("Polynomial ARD");
  setParamName("weightVariance", 0);
  weightVariance=1.0;
  addTransform(new CNegLogLogitTransform, 0);
  setParamName("biasVariance", 1);
  biasVariance=1.0;
  addTransform(new CNegLogLogitTransform, 1);
  setParamName("variance", 2);
  variance = 1.0;
  addTransform(new CNegLogLogitTransform, 2);
  scales.resize(1, getInputDim());
  gscales.resize(1, getInputDim());
  scales.setVals(0.5);
  for(int i=3; i<getInputDim()+3; i++)
    {
      addTransform(new CSigmoidTransform, i);
      string name = "inputScale";
      setParamName(name, i);
    }
  degree = 2.0;
}

double CPolyardKern::diagComputeElement(const CMatrix& X, int index) const
{
  double val=0.0;
  for(int i=0; i<getInputDim(); i++)
    {
      double x = X.getVal(index, i);
      val+=x*x*scales.getVal(i);
    }
  double arg=weightVariance*val+biasVariance;  
  return variance*pow(arg, degree);
}
// Parameters are kernel parameters
void CPolyardKern::setParam(double val, int paramNo)
{
  assert(paramNo>=0);
  assert(paramNo<nParams);
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
	  cerr << "Parameter doesn't exist.";
	  exit(1);
	}
      
    }
}
double CPolyardKern::getParam(int paramNo) const
{
  assert(paramNo>=0);
  assert(paramNo<nParams);
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
	  cerr << "Parameter doesn't exist.";
	  exit(1);
	}
    }
}
void CPolyardKern::getGradX(vector<CMatrix*>& gX, const CMatrix& X, const CMatrix& X2, const bool addG) const
{
  assert(gX.size()==X.getRows());
  assert(X.getCols()==X2.getCols());
  for(int i=0; i<X.getRows(); i++)
    {
      assert(gX[i]->getRows()==X2.getRows());
      assert(gX[i]->getCols()==X2.getCols());
      for(int k=0; k<X2.getRows(); k++)
	{
	  double valik=0.0;
	  for(int j=0; j<getInputDim(); j++)
	    {
	      valik+=X.getVal(i, j)*X2.getVal(k, j)*scales.getVal(j);
	    }
	  double arg=weightVariance*valik + biasVariance;
	  double kval=degree*variance*weightVariance*pow(arg, degree-1);
	  for(int j=0; j<X2.getCols(); j++)
	    {
	      double val = 0.0;
	      if(addG)
		val = gX[i]->getVal(k, j);
	      val += kval*X2.getVal(k, j)*scales.getVal(j);
	      gX[i]->setVal(val, k, j);
	    }
	}
    }
}
void CPolyardKern::getDiagGradX(CMatrix& gX, const CMatrix& X, const bool addG) const
{
  assert(gX.dimensionsMatch(X));
  for(int i=0; i<X.getRows(); i++)
    {
      double val=0.0;
      for(int j=0; j<getInputDim(); j++)
	{
	  double x = X.getVal(i, j);
	  val+=x*x*scales.getVal(j);
	}
      double arg=weightVariance*val + biasVariance;
      double kval=degree*variance*weightVariance*pow(arg, degree-1);
      for(int j=0; j<X.getCols(); j++)
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

double CPolyardKern::computeElement(const CMatrix& X1, int index1, 
			  const CMatrix& X2, int index2) const
{
  double valij=0.0;
  for(int i=0; i<getInputDim(); i++)
    {
      valij+=X1.getVal(index1, i)*X2.getVal(index2, i)*scales.getVal(i);
    }
  double arg=weightVariance*valij + biasVariance;
  return variance*pow(arg, degree);
}

void CPolyardKern::getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& covGrad, bool regularise) const
{
  assert(g.getRows()==1);
  assert(g.getCols()==nParams);
  double g1=0.0;
  double g2=0.0;
  double g3=0.0;
  innerProd.resize(1, X.getRows());
  gscales.zeros();
  innerProd.zeros();
  for(int i=0; i<X.getRows(); i++)
    {
      double val=0.0;
      for(int k=0; k<getInputDim(); k++)
	{
	  double x = X.getVal(i, k);
	  val+=x*x*scales.getVal(k);
	}
      innerProd.setVal(val, i);
    }
  
  for(int i=0; i<X.getRows(); i++)
    {
      for(int j=0; j<i; j++)
	{
	  double val=0.0;
	  for(int k=0; k<getInputDim(); k++)
	    {
	      val+=X.getVal(i, k)*X.getVal(j, k)*scales.getVal(k);
	    }
	  double crossProd=val;
	  double arg=weightVariance*val + biasVariance;
	  double baseCovGrad = variance*degree*pow(arg, degree-1)*covGrad.getVal(i, j);
	  g1+=baseCovGrad*val;
	  g2+=baseCovGrad;
	  g3+=pow(arg, degree)*covGrad.getVal(i, j);


	  for(int k=0; k<getInputDim(); k++)
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
  for(int i=0; i<X.getRows(); i++)
    {
      double arg = weightVariance*innerProd.getVal(i)+biasVariance;
      double baseCovGrad = variance*degree*pow(arg, degree-1)*covGrad.getVal(i, i);
      g1+=baseCovGrad*innerProd.getVal(i);
      
      g2+=baseCovGrad;
      g3+=pow(arg, degree)*covGrad.getVal(i, i);

      for(int k=0; k<getInputDim(); k++)
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
  for(int k=0; k<getInputDim(); k++)
    g.setVal(gscales.getVal(k), k+3);
  if(regularise)
    addPriorGrad(g);

}
double CPolyardKern::getGradParam(int index, const CMatrix& X, const CMatrix& covGrad) const
{
  assert(index>=0);
  assert(index<nParams);
  cerr << "Error getGradParam is not currently implemented for CPolyardKern" << endl;
  exit(1);
}

// Functions that operate on CKern.
ostream& operator<<(ostream& out, const CKern& kern)
{
  out <<  kern.display(out);
  return out;
}
void writeKernToStream(const CKern& kern, ostream& out)
{
  out << "kernVersion=" << KERNVERSION << endl;
  out << "type=" << kern.getType() << endl;
  kern.writeParamsToStream(out);
}
CKern* readKernFromStream(istream& in)
{
  CKern* pkern;
  string line;
  vector<string> tokens;
  // first line is version info.
  getline(in, line);
  ndlstrutil::tokenise(tokens, line, "=");
  if(tokens.size()!=2 || tokens[0]!="kernVersion")
    throw ndlexceptions::FileFormatError();
  if(tokens[1]!="0.1")
    throw ndlexceptions::FileVersionError();
  // next line is type.
  tokens.clear();
  getline(in, line);
  ndlstrutil::tokenise(tokens, line, "=");
  if(tokens.size()!=2 || tokens[0]!="type")
    throw ndlexceptions::FileFormatError();
  string type=tokens[1];
  if(type=="cmpnd")
    pkern = new CCmpndKern();
  else if(type=="bias")
    pkern = new CBiasKern();
  else if(type=="white")
    pkern = new CWhiteKern();
  else if(type=="lin")
    pkern = new CLinKern();
  else if(type=="linard")
    pkern = new CLinardKern();
  else if(type=="rbf")
    pkern = new CRbfKern();
  else if(type=="rbfard")
    pkern = new CRbfardKern();
  else if(type=="mlp")
    pkern = new CMlpKern();
  else if(type=="mlpard")
    pkern = new CMlpardKern();
  else if(type=="poly")
    pkern = new CPolyKern();
  else if(type=="polyard")
    pkern = new CPolyardKern();
  else
    throw ndlexceptions::FileFormatError();
  
  pkern->readParamsFromStream(in);
  return pkern;
}
void CKern::readParamsFromStream(istream& in)
{
  string line;
  vector<string> tokens;
  getline(in, line);
  ndlstrutil::tokenise(tokens, line, "=");
  if(tokens.size()>2 || tokens[0]!="numParams")
    throw ndlexceptions::FileFormatError();
  int numParams=atol(tokens[1].c_str());
  tokens.clear();
  getline(in, line);
  ndlstrutil::tokenise(tokens, line, "=");
  if(tokens.size()>2 || tokens[0]!="numFeatures")
    throw ndlexceptions::FileFormatError();
  int numFeatures=atol(tokens[1].c_str());
  CMatrix par(1, numParams);
  tokens.clear();
  getline(in, line);
  ndlstrutil::tokenise(tokens, line, " ");
  for(int i=0; i<numParams; i++)
    par.setVal(atof(tokens[i].c_str()), i);
  setInputDim(numFeatures);
  if(numParams==getNumParams())
    setParams(par);
  else
    throw ndlexceptions::FileFormatError();
  tokens.clear();
  getline(in, line);
  ndlstrutil::tokenise(tokens, line, "=");
  if(tokens.size()>2 || tokens[0]!="numPriors")
    throw ndlexceptions::FileFormatError();
  int numPriors=atol(tokens[1].c_str());
  tokens.clear();
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
  if(mxType!=type)
    {
      cerr << "Error mismatch between saved type, " << mxType << ", and Class type, " << type << "." << endl;
      exit(1);
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
void CArdKern::extractParamFromMxArray(const mxArray* matlabArray)
{
  setInputDim(mxArrayExtractIntField(matlabArray, "inputDimension"));
  assert(nParams == mxArrayExtractIntField(matlabArray, "nParams"));
  string pName;
  for(int i=0; i<nParams; i++)
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
  for(int i=0; i<nParams; i++)
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
	kernNum = addKern(new CLinKern(getInputDim()));
      else if(kernType == "rbf")  
	kernNum = addKern(new CRbfKern(getInputDim()));
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
      else
	throw ndlexceptions::Error("Unknown kernel type");
      components[kernNum]->fromMxArray(compElement);
    }
}
void CCmpndKern::addParamToMxArray(mxArray* matlabArray) const
{
  // Add comp field to mxArray.
  int dims[1];
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
#endif /* _NDLMATLAB */
