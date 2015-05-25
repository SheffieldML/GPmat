#include "CMltools.h"


CMlpMapping::CMlpMapping()
{
  _init();
}
CMlpMapping::CMlpMapping(CMatrix* pXin, CMatrix* pyOut, unsigned int nhidden, int verbos)
  :
  CMapModel(pXin->getCols(), pyOut->getCols(), pXin->getRows()), 
  COptimisable(),
  pX(pXin), py(pyOut), hiddenDim(nhidden)
{ 
  DIMENSIONMATCH(pX->getRows()==py->getRows());
  _init();
  setVerbosity(verbos);
  initStoreage();
  initVals();
}
void CMlpMapping::_init()
{
  setType("mlp");
  setName("multi-layer perceptron");
  variance = 1.0;
  setVerbosity(2);
}
void CMlpMapping::initStoreage()
{
  W1.resize(getInputDim(), hiddenDim);
  b1.resize(1, hiddenDim);
  W2.resize(hiddenDim, getOutputDim());
  b2.resize(1, getOutputDim());
  hiddenActive.resize(1, hiddenDim);
  outActive.resize(1, getOutputDim());
}
void CMlpMapping::setWeights(const CMatrix& W, unsigned int layer)
{
  BOUNDCHECK(layer>=0 && layer<2);
  switch(layer)
  {
   case 0:
     DIMENSIONMATCH(W.getRows()==getInputDim());
     DIMENSIONMATCH(W.getCols()==getHiddenDim());
     W1.deepCopy(W);
     break;
   case 1:
     DIMENSIONMATCH(W.getRows()==getHiddenDim());
     DIMENSIONMATCH(W.getCols()==getOutputDim());
     W2.deepCopy(W);
     break;
   default:
     throw ndlexceptions::Error("Unknown layer");
  }
}
void CMlpMapping::setBias(const CMatrix& b, unsigned int layer)
{
  BOUNDCHECK(layer>=0 && layer<2);
  switch(layer)
  {
  case 0:
    DIMENSIONMATCH(b.getRows() == 1);
    DIMENSIONMATCH(b.getCols() == getHiddenDim());
    b1.deepCopy(b);
    break;
  case 1:
    DIMENSIONMATCH(b.getRows() == 1);
    DIMENSIONMATCH(b.getCols() == getOutputDim());
    b2.deepCopy(b);
    break;
  default:
     throw ndlexceptions::Error("Unknown layer");
  }
}
void CMlpMapping::initVals() 
{
  W1.randn(1.0/(double)(getInputDim()+1), 0.0);
  b1.randn(1.0/(double)(getInputDim()+1), 0.0);
  W2.randn(1.0/(double)(hiddenDim+1), 0.0);
  b2.randn(1.0/(double)(hiddenDim+1), 0.0);
}


unsigned int CMlpMapping::getOptNumParams() const
{
  return (getInputDim()+1)*hiddenDim + (hiddenDim+1)*getOutputDim();
}

void CMlpMapping::getOptParams(CMatrix& param) const
{
  int counter = 0;
  for(unsigned int j=0; j<hiddenDim; j++)
  {
    for(unsigned int i=0; i<getInputDim(); i++)
    {
      param.setVal(W1.getVal(i, j), counter);
      counter++;
    }
  }
  for(unsigned int j=0; j<hiddenDim; j++)
  {
    param.setVal(b1.getVal(j), counter);
    counter++;
  }
  for(unsigned int j=0; j<getOutputDim(); j++)
  {
    for(unsigned int i=0; i<hiddenDim; i++)
    {
      param.setVal(W2.getVal(i, j), counter);
      counter++;
    }
  }
  for(unsigned int j=0; j<getOutputDim(); j++)
  {
    param.setVal(b2.getVal(j), counter);
    counter++;
  }
}

void CMlpMapping::setOptParams(const CMatrix& param)
{
  int counter = 0;
  for(unsigned int j=0; j<hiddenDim; j++)
  {
    for(unsigned int i=0; i<getInputDim(); i++)
    {
      W1.setVal(param.getVal(counter), i, j);
      counter++;
    }
  }
  for(unsigned int j=0; j<hiddenDim; j++)
  {
    b1.setVal(param.getVal(counter), j);
    counter++;
  }
  for(unsigned int j=0; j<getOutputDim(); j++)
  {
    for(unsigned int i=0; i<hiddenDim; i++)
    {
      W2.setVal(param.getVal(counter), i, j);
      counter++;
    }
  }
  for(unsigned int j=0; j<getOutputDim(); j++)
  {
    b2.setVal(param.getVal(counter), j);
    counter++;
  }
}
double CMlpMapping::outGradParams(CMatrix& g, const CMatrix& Xin, unsigned int pointNo, unsigned int outputNo) const
{

  hiddenActive.deepCopy(b1);
  hiddenActive.gemvRowRow(0, W1, Xin, pointNo, 1.0, 1.0, "t");
  hiddenActive.tanh();
  tanhActive.deepCopy(hiddenActive);
  tanhActive.dotMultiply(tanhActive);
  tanhActive.negate();
  tanhActive.add(1.0);
  int counter = 0;
  for(unsigned int j=0; j<hiddenDim; j++)
  {
    for(unsigned int i=0; i<getInputDim(); i++)
    {
      g.setVal(Xin.getVal(pointNo, i)*tanhActive.getVal(j)*W2.getVal(j, outputNo), counter);
      counter++;
    }
  }
  for(unsigned int j=0; j<hiddenDim; j++)
  {
    g.setVal(tanhActive.getVal(j)*W2.getVal(j, outputNo), counter);
    counter++;
  }
  for(unsigned int j=0; j<getOutputDim(); j++)
  { 
    if(j==outputNo)
    {
      for(unsigned int i=0; i<hiddenDim; i++)
      {
        g.setVal(hiddenActive.getVal(i), counter);
        counter++;
      }
    }
    else
    {
      for(unsigned int i=0; i<hiddenDim; i++)
      {
        g.setVal(0.0, counter);
        counter++;
      }
    }
  }
  for(unsigned int j=0; j<getOutputDim(); j++)
  {
    if(j==outputNo)
    {
      g.setVal(1.0, counter);
      counter++;
    }
    else
    {
      g.setVal(0.0, counter);
    }
  }
 double out = b2.getVal(outputNo);
 out += hiddenActive.dotRowCol(0, W2, outputNo);
 return out;

}
void CMlpMapping::out(CMatrix& yPred, const CMatrix& Xin) const
{
  DIMENSIONMATCH(yPred.getRows()==Xin.getRows());
  DIMENSIONMATCH(Xin.getCols()==getInputDim());
  for(unsigned int i=0; i<yPred.getRows(); i++)
  {
    hiddenActive.deepCopy(b1);
    hiddenActive.gemvRowRow(0, W1, Xin, i, 1.0, 1.0, "t");
    hiddenActive.tanh();
    yPred.copyRowRow(i, b2, 0);
    yPred.gemvRowRow(i, W2, hiddenActive, 0, 1.0, 1.0, "t");
  }

}
double CMlpMapping::outGradX(CMatrix& g, const CMatrix &Xin, unsigned int pointNo, unsigned int outputNo) const
{
  throw ndlexceptions::NotImplementedError("outGradX not yet implemented for CMlpMapping.");
}

// compute the log likelihood.
double CMlpMapping::logLikelihood() const
{
  double L=0.0;
  for(unsigned int i=0; i<getNumData(); i++)
  {
    hiddenActive.deepCopy(b1);
    hiddenActive.gemvRowRow(0, W1, *pX, i, 1.0, 1.0, "t");
    hiddenActive.tanh();
    outActive.copyRowRow(0, b2, 0);
    outActive.gemvRowRow(0, W2, hiddenActive, 0, 1.0, 1.0, "t");
    L+=outActive.dist2Row(0, *py, i);
  }
  L=L/variance;
  L+=(double)getNumData()*(ndlutil::LOGTWOPI + log(variance)); 
  L*=-0.5;
  //L+=priorLogProb();
  return L;
}
// compute the gradients wrt parameters and latent variables.
double CMlpMapping::logLikelihoodGradient(CMatrix& g) const 
{
  double L=0.0;
  g.zeros();
  CMatrix gtemp(1, getOptNumParams());
  for(unsigned int j=0; j<getOutputDim(); j++)
  {
    for(unsigned int i=0; i<getNumData(); i++)
    {
      double diff = outGradParams(gtemp, *pX, i, j) - py->getVal(i, j);
      L+= diff*diff;
      gtemp.scale(-diff/variance);
      g.add(gtemp);
    }
  }

  L=L/variance;
  L+=(double)getNumData()*(ndlutil::LOGTWOPI + log(variance)); 
  L*=-0.5;
  //L+=priorLogProb();
  return L;
}
void CMlpMapping::updateG() const 
{
}

double CMlpMapping::pointLogLikelihood(const CMatrix& yTest, const CMatrix& Xtest) const
{
  out(outActive, Xtest);
  double L = outActive.dist2Row(0, yTest, 0);
  L = L/variance;
  L += ndlutil::LOGTWOPI + log(variance);
  L*=-0.5;
  return L;
}
#ifdef _NDLMATLAB
CMlpMapping::CMlpMapping(CMatrix* inData, 
	   CMatrix* targetData, 
	   const string mlpMappingInfoFile, 
	   const string mlpMappingInfoVariable, 
	   int verbos) : 
  CMapModel(inData->getCols(), targetData->getCols(), inData->getRows()), COptimisable(),
  pX(inData), py(targetData)
{
  _init();
  setVerbosity(verbos);
  readMatlabFile(mlpMappingInfoFile, mlpMappingInfoVariable);
}
mxArray* CMlpMapping::toMxArray() const
{
  int dims[1];
  dims[0]=1;
  mxArray* matlabArray;
  const char* fieldNames[]={"type", "nin", "nhidden", "nout", "nwts", "outfn", "w1", "b1", "w2", "b2", "numParams", "hiddenDim", "inputDim", "outputDim", "optimiser"};
  matlabArray = mxCreateStructArray(1, dims, 15, fieldNames);
  mxSetField(matlabArray, 0, "type", convertMxArray("mlp"));
  mxSetField(matlabArray, 0, "nin", convertMxArray((double)getInputDim()));
  mxSetField(matlabArray, 0, "nhidden", convertMxArray((double)hiddenDim));
  mxSetField(matlabArray, 0, "nout", convertMxArray((double)getOutputDim()));
  mxSetField(matlabArray, 0, "nwts", convertMxArray((double)getOptNumParams()));
  mxSetField(matlabArray, 0, "outfn", convertMxArray("linear"));
  mxSetField(matlabArray, 0, "w1", convertMxArray("W1"));
  mxSetField(matlabArray, 0, "b1", convertMxArray("b1"));
  mxSetField(matlabArray, 0, "w2", convertMxArray("W2"));
  mxSetField(matlabArray, 0, "b2", convertMxArray("b2"));
  mxSetField(matlabArray, 0, "numParams", convertMxArray((double)getOptNumParams()));
  mxSetField(matlabArray, 0, "hiddenDim", convertMxArray((double)hiddenDim));
  mxSetField(matlabArray, 0, "inputDim", convertMxArray((double)getInputDim()));
  mxSetField(matlabArray, 0, "outputDim", convertMxArray((double)getOutputDim()));
  mxSetField(matlabArray, 0, "optimiser", convertMxArray(getDefaultOptimiserStr()));
  return matlabArray;
}
void CMlpMapping::fromMxArray(const mxArray* matlabArray)
{
  string mxType = mxArrayExtractStringField(matlabArray, "type");
  if(mxType!=getType())
  {
    throw ndlexceptions::FileReadError("Error mismatch between saved type, " + mxType + ", and Class type, " + getType() + ".");
  }
  setInputDim(mxArrayExtractIntField(matlabArray, "inputDim"));
  setOutputDim(mxArrayExtractIntField(matlabArray, "outputDim"));
  hiddenDim = mxArrayExtractIntField(matlabArray, "hiddenDim");
  setDefaultOptimiserStr(mxArrayExtractStringField(matlabArray, "optimiser"));
  initStoreage();
  W1.fromMxArray(mxArrayExtractMxArrayField(matlabArray, "w1"));
  b1.fromMxArray(mxArrayExtractMxArrayField(matlabArray, "b1"));
  W2.fromMxArray(mxArrayExtractMxArrayField(matlabArray, "w2"));
  b2.fromMxArray(mxArrayExtractMxArrayField(matlabArray, "b2"));
}

#else /* not _NDLMATLAB */
#endif
// Optimise the mapping with respect to the parameters.
void CMlpMapping::optimise(unsigned int iters)
{
  if(getVerbosity()>2)
  {
    cout << "Initial model:" << endl;
    display(cout);
  }
  if(getVerbosity()>2 && getOptNumParams()<40)
    checkGradients();
  setMaxIters(iters);
  runDefaultOptimiser();
  
  if(getVerbosity()>1)
    cout << "... done. " << endl;
  if(getVerbosity()>0)
    display(cout);
}

bool CMlpMapping::equals(const CMlpMapping& model, const double tol) const 
{
  if(model.getHiddenDim()!=getHiddenDim())
    return false;
  if(model.getInputDim()!=getInputDim())
    return false;
  if(model.getOutputDim()!=getOutputDim())
    return false;
  if(!W1.equals(model.W1, tol))
    return false;
  if(!b1.equals(model.b1, tol))
    return false;
  if(!W2.equals(model.W2, tol))
    return false;
  if(!b2.equals(model.b2, tol))
    return false;
  return true;
}

void CMlpMapping::display(ostream& os) const 
{
  cout << "Multi-Layer Perceptron Model:" << endl;
  cout << "Optimiser: " << getDefaultOptimiserStr() << endl;
  cout << "Data Set Size: " << getNumData() << endl;
  cout << "Number hidden: " << hiddenDim << endl;
  cout << "Log likelihood: " << logLikelihood() << endl;
}


void CMlpMapping::writeParamsToStream(ostream& out) const
{
  writeToStream(out, "baseType", getBaseType());
  writeToStream(out, "type", getType());
  writeToStream(out, "numData", getNumData());
  writeToStream(out, "outputDim", getOutputDim());
  writeToStream(out, "inputDim", getInputDim());

  writeToStream(out, "hiddenDim", hiddenDim);
  writeToStream(out, "numParams", getOptNumParams());
  CMatrix par(1, getOptNumParams());
  getOptParams(par);
  par.toStream(out);
}

void CMlpMapping::readParamsFromStream(istream& in)
{
  string tbaseType = getBaseTypeStream(in);
  if(tbaseType != getBaseType())
    throw ndlexceptions::StreamFormatError("baseType", "Error mismatch between saved base type, " + tbaseType + ", and Class base type, " + getType() + ".");
  string ttype = getTypeStream(in);
  if(ttype != getType())
    throw ndlexceptions::StreamFormatError("type", "Error mismatch between saved type, " + ttype + ", and Class type, " + getType() + ".");
  setNumData(readIntFromStream(in, "numData"));
  setOutputDim(readIntFromStream(in, "outputDim"));
  setInputDim(readIntFromStream(in, "inputDim"));

  hiddenDim = readIntFromStream(in, "hiddenDim");
  initStoreage();

  unsigned int nPar = readIntFromStream(in, "numParams");
  if(nPar!=getOptNumParams())
    throw ndlexceptions::StreamFormatError("numParams", "Number of parameters does not match.");
  CMatrix par(1, getOptNumParams());
  par.fromStream(in);
  setOptParams(par);
}

// Functions which operate on the object
void writeMlpMappingToStream(const CMlpMapping& model, ostream& out) 
{
  model.toStream(out);
}

void writeMlpMappingToFile(const CMlpMapping& model, const string modelFileName, const string comment)
{
  if(model.getVerbosity()>0) {
    cout << "Saving model file." << endl;
  }
  ofstream out(modelFileName.c_str());
  if(!out) throw ndlexceptions::FileWriteError(modelFileName);
  if(comment.size()>0)
    out << "# " << comment << endl;
  writeMlpMappingToStream(model, out);
  out.close();

}

CMlpMapping* readMlpMappingFromStream(istream& in)
{
  CMlpMapping* pmodel = new CMlpMapping();
  pmodel->fromStream(in);
  return pmodel;
}
    
CMlpMapping* readMlpMappingFromFile(const string modelFileName, int verbosity)
{
  // File is m, beta, X
  if(verbosity>0)
    cout << "Loading model file." << endl;
  ifstream in(modelFileName.c_str());
  if(!in.is_open()) throw ndlexceptions::FileReadError(modelFileName);
  CMlpMapping* pmodel;
  try 
  {
    pmodel = readMlpMappingFromStream(in);
  }
  catch(ndlexceptions::StreamFormatError err) 
  {
    throw ndlexceptions::FileFormatError(modelFileName, err);
  }
  if(verbosity>0)
    cout << "... done." << endl;
  in.close();
  return pmodel;

}

CLinearMapping::CLinearMapping()
{
  _init();
}
CLinearMapping::CLinearMapping(CMatrix* pXin, CMatrix* pyOut, int verbos)
  :
  CMapModel(pXin->getCols(), pyOut->getCols(), pXin->getRows()), 
  COptimisable(),
  pX(pXin), py(pyOut)
{ 
  DIMENSIONMATCH(pX->getRows()==py->getRows());
  _init();
  setVerbosity(verbos);
  initStoreage();
  initVals();
  
}
void CLinearMapping::_init()
{
  setType("linear");
  setName("linear mapping");
  variance = 1.0;
  setVerbosity(2);
}
void CLinearMapping::initStoreage()
{
  W.resize(getInputDim(), getOutputDim());
  b.resize(1, getOutputDim());
  outActive.resize(1, getOutputDim());
}
void CLinearMapping::setWeights(const CMatrix& Win, unsigned int layer)
{
  BOUNDCHECK(layer==0);
  switch(layer)
  {
   case 0:
     DIMENSIONMATCH(Win.getRows()==getInputDim());
     DIMENSIONMATCH(Win.getCols()==getOutputDim());
     W.deepCopy(Win);
     break;
   default:
     throw ndlexceptions::Error("Unknown layer");
  }
}
void CLinearMapping::setBias(const CMatrix& bin, unsigned int layer)
{
  BOUNDCHECK(layer>=0 && layer<2);
  switch(layer)
  {
  case 0:
    DIMENSIONMATCH(bin.getRows() == 1);
    DIMENSIONMATCH(bin.getCols() == getOutputDim());
    b.deepCopy(bin);
    break;
  default:
     throw ndlexceptions::Error("Unknown layer");
  }
}
void CLinearMapping::initVals() 
{
  W.randn(1.0/(double)(getInputDim()+1), 0.0);
  b.randn(1.0/(double)(getInputDim()+1), 0.0);
}


unsigned int CLinearMapping::getOptNumParams() const
{
  return (getInputDim()+1)*getOutputDim();
}

void CLinearMapping::getOptParams(CMatrix& param) const
{
  int counter = 0;
  for(unsigned int j=0; j<getOutputDim(); j++)
  {
    for(unsigned int i=0; i<getInputDim(); i++)
    {
      param.setVal(W.getVal(i, j), counter);
      counter++;
    }
  }
  for(unsigned int j=0; j<getOutputDim(); j++)
  {
    param.setVal(b.getVal(j), counter);
    counter++;
  }
}

void CLinearMapping::setOptParams(const CMatrix& param)
{
  int counter = 0;
  for(unsigned int j=0; j<getOutputDim(); j++)
  {
    for(unsigned int i=0; i<getInputDim(); i++)
    {
      W.setVal(param.getVal(counter), i, j);
      counter++;
    }
  }
  for(unsigned int j=0; j<getOutputDim(); j++)
  {
    b.setVal(param.getVal(counter), j);
    counter++;
  }
}
double CLinearMapping::outGradParams(CMatrix& g, const CMatrix& Xin, unsigned int pointNo, unsigned int outputNo) const
{

  int counter = 0;
  for(unsigned int j=0; j<getOutputDim(); j++)
  { 
    if(j==outputNo)
    {
      for(unsigned int i=0; i<getInputDim(); i++)
      {
        g.setVal(Xin.getVal(pointNo, i), counter);
        counter++;
      }
    }
    else
    {
      for(unsigned int i=0; i<getInputDim(); i++)
      {
        g.setVal(0.0, counter);
        counter++;
      }
    }
  }
  for(unsigned int j=0; j<getOutputDim(); j++)
  {
    if(j==outputNo)
    {
      g.setVal(1.0, counter);
      counter++;
    }
    else
    {
      g.setVal(0.0, counter);
    }
  }
 double out = b.getVal(outputNo);
 out += Xin.dotRowCol(pointNo, W, outputNo);
 return out;

}
void CLinearMapping::out(CMatrix& yPred, const CMatrix& Xin) const
{
  DIMENSIONMATCH(yPred.getRows()==Xin.getRows());
  DIMENSIONMATCH(Xin.getCols()==getInputDim());
  for(unsigned int i=0; i<yPred.getRows(); i++)
  {
    yPred.copyRowRow(i, b, 0);
    yPred.gemvRowRow(i, W, Xin, i, 1.0, 1.0, "t");
  }

}
double CLinearMapping::outGradX(CMatrix& g, const CMatrix &Xin, unsigned int pointNo, unsigned int outputNo) const
{
  throw ndlexceptions::NotImplementedError("outGradX not yet implemented for CLinearMapping.");
}

// compute the log likelihood.
double CLinearMapping::logLikelihood() const
{
  double L=0.0;
  for(unsigned int i=0; i<getNumData(); i++)
  {
    outActive.copyRowRow(0, b, 0);
    outActive.gemvRowRow(0, W, *pX, i, 1.0, 1.0, "t");
    L+=outActive.dist2Row(0, *py, i);
  }
  L=L/variance;
  L+=(double)getNumData()*(ndlutil::LOGTWOPI + log(variance)); 
  L*=-0.5;
  //L+=priorLogProb();
  return L;
}
// compute the gradients wrt parameters and latent variables.
double CLinearMapping::logLikelihoodGradient(CMatrix& g) const 
{
  double L=0.0;
  g.zeros();
  CMatrix gtemp(1, getOptNumParams());
  for(unsigned int j=0; j<getOutputDim(); j++)
  {
    for(unsigned int i=0; i<getNumData(); i++)
    {
      double diff = outGradParams(gtemp, *pX, i, j) - py->getVal(i, j);
      L+= diff*diff;
      gtemp.scale(-diff/variance);
      g.add(gtemp);
    }
  }

  L=L/variance;
  L+=(double)getNumData()*(ndlutil::LOGTWOPI + log(variance)); 
  L*=-0.5;
  //L+=priorLogProb();
  return L;
}
void CLinearMapping::updateG() const 
{
}

double CLinearMapping::pointLogLikelihood(const CMatrix& yTest, const CMatrix& Xtest) const
{
  out(outActive, Xtest);
  double L = outActive.dist2Row(0, yTest, 0);
  L = L/variance;
  L += ndlutil::LOGTWOPI + log(variance);
  L*=-0.5;
  return L;
}
#ifdef _NDLMATLAB
CLinearMapping::CLinearMapping(CMatrix* inData, 
	   CMatrix* targetData, 
	   const string linearMappingInfoFile, 
	   const string linearMappingInfoVariable, 
	   int verbos) : 
  CMapModel(inData->getCols(), targetData->getCols(), inData->getRows()), COptimisable(),
  pX(inData), py(targetData)
{
  _init();
  setVerbosity(verbos);
  readMatlabFile(linearMappingInfoFile, linearMappingInfoVariable);
}
mxArray* CLinearMapping::toMxArray() const
{
  int dims[1];
  dims[0]=1;
  mxArray* matlabArray;
  const char* fieldNames[]={"type", "W", "b", "numParams", "inputDim", "outputDim", "optimiser", "beta"};
  matlabArray = mxCreateStructArray(1, dims, 8, fieldNames);
  mxSetField(matlabArray, 0, "type", convertMxArray("linear"));
  mxSetField(matlabArray, 0, "W", convertMxArray("W"));
  mxSetField(matlabArray, 0, "b", convertMxArray("b"));
  mxSetField(matlabArray, 0, "numParams", convertMxArray((double)getOptNumParams()));
  mxSetField(matlabArray, 0, "inputDim", convertMxArray((double)getInputDim()));
  mxSetField(matlabArray, 0, "outputDim", convertMxArray((double)getOutputDim()));
  mxSetField(matlabArray, 0, "optimiser", convertMxArray(getDefaultOptimiserStr()));
  mxSetField(matlabArray, 0, "beta", convertMxArray(variance));
  return matlabArray;
}
void CLinearMapping::fromMxArray(const mxArray* matlabArray)
{
  string mxType = mxArrayExtractStringField(matlabArray, "type");
  if(mxType!=getType())
  {
    throw ndlexceptions::FileReadError("Error mismatch between saved type, " + mxType + ", and Class type, " + getType() + ".");
  }
  setInputDim(mxArrayExtractIntField(matlabArray, "inputDim"));
  setOutputDim(mxArrayExtractIntField(matlabArray, "outputDim"));
  //setDefaultOptimiserStr(mxArrayExtractStringField(matlabArray, "optimiser"));
  initStoreage();
  W.fromMxArray(mxArrayExtractMxArrayField(matlabArray, "W"));
  b.fromMxArray(mxArrayExtractMxArrayField(matlabArray, "b"));
}

#else /* not _NDLMATLAB */
#endif
// Optimise the mapping
void CLinearMapping::optimise(unsigned int iters)
{
  if(getVerbosity()>2)
  {
    cout << "Initial model:" << endl;
    display(cout);
  }
  if(getVerbosity()>2 && getOptNumParams()<40)
    checkGradients();
  setMaxIters(iters);
  runDefaultOptimiser();
  
  if(getVerbosity()>1)
    cout << "... done. " << endl;
  if(getVerbosity()>0)
    display(cout);
}

bool CLinearMapping::equals(const CLinearMapping& model, const double tol) const 
{
  if(model.getInputDim()!=getInputDim())
    return false;
  if(model.getOutputDim()!=getOutputDim())
    return false;
  if(!W.equals(model.W, tol))
    return false;
  if(!b.equals(model.b, tol))
    return false;
  return true;
}

void CLinearMapping::display(ostream& os) const 
{
  cout << "Linear Mapping:" << endl;
  cout << "Optimiser: " << getDefaultOptimiserStr() << endl;
  cout << "Data Set Size: " << getNumData() << endl;
  cout << "Log likelihood: " << logLikelihood() << endl;
}


void CLinearMapping::writeParamsToStream(ostream& out) const
{
  writeToStream(out, "baseType", getBaseType());
  writeToStream(out, "type", getType());
  writeToStream(out, "numData", getNumData());
  writeToStream(out, "outputDim", getOutputDim());
  writeToStream(out, "inputDim", getInputDim());

  writeToStream(out, "numParams", getOptNumParams());
  CMatrix par(1, getOptNumParams());
  getOptParams(par);
  par.toStream(out);
}

void CLinearMapping::readParamsFromStream(istream& in)
{
  string tbaseType = getBaseTypeStream(in);
  if(tbaseType != getBaseType())
    throw ndlexceptions::StreamFormatError("baseType", "Error mismatch between saved base type, " + tbaseType + ", and Class base type, " + getType() + ".");
  string ttype = getTypeStream(in);
  if(ttype != getType())
    throw ndlexceptions::StreamFormatError("type", "Error mismatch between saved type, " + ttype + ", and Class type, " + getType() + ".");
  setNumData(readIntFromStream(in, "numData"));
  setOutputDim(readIntFromStream(in, "outputDim"));
  setInputDim(readIntFromStream(in, "inputDim"));
  initStoreage();

  unsigned int nPar = readIntFromStream(in, "numParams");
  if(nPar!=getOptNumParams())
    throw ndlexceptions::StreamFormatError("numParams", "Number of parameters does not match.");
  CMatrix par(1, getOptNumParams());
  par.fromStream(in);
  setOptParams(par);
}

// Functions which operate on the object
void writeLinearMappingToStream(const CLinearMapping& model, ostream& out) 
{
  model.toStream(out);
}

void writeLinearMappingToFile(const CLinearMapping& model, const string modelFileName, const string comment)
{
  if(model.getVerbosity()>0) {
    cout << "Saving model file." << endl;
  }
  ofstream out(modelFileName.c_str());
  if(!out) throw ndlexceptions::FileWriteError(modelFileName);
  if(comment.size()>0)
    out << "# " << comment << endl;
  writeLinearMappingToStream(model, out);
  out.close();

}

CLinearMapping* readLinearMappingFromStream(istream& in)
{
  CLinearMapping* pmodel = new CLinearMapping();
  pmodel->fromStream(in);
  return pmodel;
}
    
CLinearMapping* readLinearMappingFromFile(const string modelFileName, int verbosity)
{
  // File is m, beta, X
  if(verbosity>0)
    cout << "Loading model file." << endl;
  ifstream in(modelFileName.c_str());
  if(!in.is_open()) throw ndlexceptions::FileReadError(modelFileName);
  CLinearMapping* pmodel;
  try 
  {
    pmodel = readLinearMappingFromStream(in);
  }
  catch(ndlexceptions::StreamFormatError err) 
  {
    throw ndlexceptions::FileFormatError(modelFileName, err);
  }
  if(verbosity>0)
    cout << "... done." << endl;
  in.close();
  return pmodel;

}

