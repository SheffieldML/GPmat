#include "CIvm.h"

CIvm::CIvm(const CMatrix& inData, const CMatrix& targetData, 
	   CKern& kernel, CNoise& noiseModel, const int selectCrit,
	   const int dVal, const int verbos) 
  : X(inData), y(targetData), kern(kernel), 
    noise(noiseModel), selectionCriterion(selectCrit), numTarget(y.getCols()), numData(y.getRows()), lastEntropyChange(0.0), cumEntropy(0.0), activeSetSize(dVal), verbosity(verbos)
{
  assert(dVal<=numData);
  assert(X.getRows()==y.getRows());
  terminate = false;
  noise.setVerbosity(getVerbosity());
  epUpdate = false;
  // need to set up Sigma for handing spherical and non-spherical.
  // set numCovStruct
  init();
}

CIvm::CIvm(const CMatrix& inData, 
	   const CMatrix& targetData, 
	   CKern& kernel, 
	   CNoise& noiseModel, 
	   const string ivmInfoFile, 
	   const string ivmInfoVariable, 
	   const int verbos) : 
  X(inData), y(targetData), 
  kern(kernel), noise(noiseModel), 
  numTarget(y.getCols()), numData(y.getRows()),  
  lastEntropyChange(0.0), cumEntropy(0.0), 
  verbosity(verbos), epUpdate(false), terminate(false)
{
  readMatlabFile(ivmInfoFile, ivmInfoVariable);
  initStoreage(); // storeage has to be allocated after finding active set size.
  for(int i=0; i<numData; i++)
    for(int j=0; j<activeSetSize; j++)
      Kstore.setVal(kern.computeElement(X, i, X, activeSet[j]), i, j);
  
  Kstore.getMatrix(K, activeSet, 0, activeSetSize-1);

  for(int i=0; i<activeSet.size(); i++)
    K.setVal(kern.diagComputeElement(X, activeSet[i]), i, i);
  K.writeMatlabFile("crap.mat", "K");
  for(int j=0; j<numCovStruct; j++)
    {
      L[j].deepCopy(K);
      L[j].updateMatlabFile("crap.mat", "KL");
      for(int i=0; i<activeSetSize; i++)
	{
	  double lval = L[j].getVal(i, i);
	  lval += 1/beta.getVal(i, j);
	  L[j].setVal(lval, i, i);
	}
      L[j].updateMatlabFile("crap.mat", "KB");
      L[j].setSymmetric(true);
      L[j].chol("L");
      L[j].updateMatlabFile("crap.mat", "L");
      Linv[j].deepCopy(L[j]);
      // TODO should not use regular inverse here as matrix is lower triangular.
      Linv[j].inv();
      M[j].gemm(Linv[j], Kstore, 1.0, 0.0, "n", "t");
    }
  for(int i=0; i<activeSetSize; i++)
    activeX.copyRowRow(i, X, activeSet[i]);
  CMatrix mu(numData, numTarget);
  CMatrix varSigma(numData, numTarget);
  posteriorMeanVar(mu, varSigma, X);
  noise.setMus(mu);
  noise.setVarSigmas(varSigma);
  updateNuG();
}
void CIvm::posteriorMeanVar(CMatrix& mu, CMatrix& varSigma, const CMatrix& X) const
{
  assert(mu.getCols()==numTarget);
  assert(varSigma.getCols()==numTarget);
  CMatrix kX(X.getRows(), activeSetSize);
  kern.compute(kX, X, activeX);
  if(numCovStruct==1)
    {
      kX.trsm(L[0], 1.0, "R", "L", "T", "N"); // now it is Linvk
      for(int i=0; i<X.getRows(); i++)
	{
	  double vsVal = kern.diagComputeElement(X, i) - kX.norm2Row(i);
	  assert(vsVal>=0);
	  for(int j=0; j<numTarget; j++)
	    varSigma.setVal(vsVal, i, j);
	  
	}
      kX.trsm(L[0], 1.0, "R", "L", "N", "N"); // now it is Kinvk
      mu.gemm(kX, m, 1.0, 0.0, "N", "N");
    }
  else
    {
      CMatrix Lk(X.getRows(), activeSetSize);
      
      for(int k=0; k<numCovStruct; k++)
	{
	  Lk.deepCopy(kX);
	  Lk.trsm(L[k], 1.0, "R", "L", "T", "N");
	  for(int i=0; i<X.getRows(); i++)
	    {
	      double vsVal=kern.diagComputeElement(X, i) - Lk.norm2Row(i);
	      assert(vsVal>=0);
	      varSigma.setVal(vsVal, i, k);
	    }
	  Lk.trsm(L[k], 1.0, "R", "L", "N", "N"); // now it is Kinvk
	  mu.gemvColCol(k, Lk, m, k, 1.0, 0.0, "N");
	}
    }
}
void CIvm::init()
{
  initStoreage();
  initVals();
}
void CIvm::initStoreage()
{
  if(noise.isSpherical())
    numCovStruct = 1;
  else
    numCovStruct = numTarget;
  Kstore.resize(numData, activeSetSize);
  m.resize(activeSetSize, y.getCols());
  beta.resize(activeSetSize, y.getCols());
  nu.resize(y.getRows(), y.getCols());
  g.resize(y.getRows(), y.getCols());

  // set up s, a and ainv
  s.resize(numData, 1); // s is a columnrow vector.
  a.resize(activeSetSize, 1); // a is a column vector.
  ainv.resize(1, activeSetSize); // ainv is a row vector.
  
  // set up L, M and Linv
  M = new CMatrix[numCovStruct];
  L = new CMatrix[numCovStruct];
  Linv = new CMatrix[numCovStruct];
  for(int c=0; c<numCovStruct; c++)
    {
      M[c].resize(activeSetSize, numData);
      L[c].resize(activeSetSize, activeSetSize);
      Linv[c].resize(activeSetSize, activeSetSize);
    }
  // set up K invK and covGrad
  K.resize(activeSetSize, activeSetSize);
  activeX.resize(activeSetSize, X.getCols());
  invK.resize(activeSetSize, activeSetSize);
  invK.setSymmetric(true);
  covGrad.resize(activeSetSize, activeSetSize);
  covGrad.setSymmetric(true);
}

void CIvm::initVals()
{
  // set Kstore to zeros numData, activeSetSize.
  Kstore.setVals(0.0);
  // set m and beta to zeros(size of y) -- could do with sparse representation.
  m.setVals(0.0);
  beta.setVals(0.0);
  // set nu to zeros(size of y)
  nu.setVals(0.0);
  // set g to zeros(size of y)
  g.setVals(0.0);
  // set noise.varSigma to diagonal of kernel.
  noise.setMus(0.0);
  double dk=0.0;
  for(int i=0; i<y.getRows(); i++)
    {
      dk = kern.diagComputeElement(X, i);
      for(int j=0; j<y.getCols(); j++)
	{
	  noise.setVarSigma(dk, i, j);
	}
    }
  for(int c=0; c<numCovStruct; c++)
    {
      M[c].setVals(0.0);
      L[c].setVals(0.0);
      Linv[c].setVals(0.0);
    }
  invK.zeros();
  covGrad.zeros();
  
  // fill the inactive set.
  inactiveSet.erase(inactiveSet.begin(), inactiveSet.end());
  for(int i=0; i<numData; i++)
    {
      inactiveSet.push_back(i);
    }
  
  // empty the active set.
  activeSet.erase(activeSet.begin(), activeSet.end());
  updateNuG();
}
void CIvm::selectPoints()
{
  int index=0;
  if(getVerbosity()>1)
    cout << "Selecting " << activeSetSize << " points ... " << endl;
  for(int k=0; k<activeSetSize; k++)
    {
      index = selectPointAdd();
      addPoint(index);
      if(getVerbosity()>2)
	cout << k << "th addition: added point " << index << endl;
      
    }
  if(getVerbosity()>1)
    cout << "... done." << endl;
  if(isEpUpdate())
    cerr << "EP update not yet implemented.";
}
void CIvm::addPoint(const int index)
{
  // check index is in inactive set
  assert(find(inactiveSet.begin(), inactiveSet.end(), index)!=inactiveSet.end());
  vector<int>::iterator pos = find(inactiveSet.begin(), inactiveSet.end(), index);
  updateSite(index);
  updateM(index);
  inactiveSet.erase(pos);
  activeX.copyRowRow(activeSet.size(), X, index);
  activeSet.push_back(index);
  updateNuG();
}
void CIvm::updateSite(const int index)
{
  int actIndex = activeSet.size();
  noise.updateSites(m, beta, actIndex, g, nu, index);
  for(int j=0; j<beta.getCols(); j++)
    {
      double betVal = beta.getVal(actIndex, j);
      if(betVal<0)
	{
	  if(noise.isLogConcave())
	    {
	      cerr << "Error beta less than zero for log concave model.";
	    }
	  else
	    {
	      beta.setVal(1e-6, actIndex, j);
	      cout << "Beta less than zero fixing to 1e-6." << endl;
	    }
	}
    }
  
}

void CIvm::updateM(const int index)
{
  int activePoint = activeSet.size();
  for(int i=0; i<Kstore.getRows(); i++)
    {
      Kstore.setVal(kern.computeElement(X, i, X, index), i, activePoint);      
    }
  // add white noise term to relevant index.
  double val = Kstore.getVal(index, activePoint);
  Kstore.setVal(val+kern.getWhite(), index, activePoint);
  double lValInv = 0.0;
  double vs = 0.0;
  double ms = 0.0;
  double sVal = 0.0;
  for(int c=0; c<numCovStruct; c++)
    {      
      lValInv = sqrt(nu.getVal(index, c));
      // set s from the kernel -- it is a column vector..
      Kstore.getMatrix(s, 0, numData-1, activePoint, activePoint);
      M[c].getMatrix(a, 0, activeSetSize-1, index, index);
      s.gemv(M[c], a, -1.0, 1.0, "t");
      if(lValInv<NULOW)
	cout << "Warning: square root of nu is " << lValInv << endl;
      // place the vector s at the bottom of M.
      s.trans();
      M[c].setMatrix(activePoint, 0, s);
      s.trans();
      M[c].scaleRow(activePoint, lValInv);
      a.trans(); // turn a into a row vector.
      L[c].setMatrix(activePoint, 0, a);
      L[c].setVal(1/lValInv, activePoint, activePoint);
      a.trans(); // turn a into a column vector.
      // update the varSigma and mu systems.
      double varSig = 0.0;
      for(int i=0; i<numData; i++)
	{
	  sVal = s.getVal(i, 0);
	  varSig = noise.getVarSigma(i, c)
	    -sVal*sVal*nu.getVal(index, c);
	  if(isnan(varSig))
	    cout << "varSigma is varSig" << endl;
	  noise.setVarSigma(varSig, i, c);
	  noise.setMu(noise.getMu(i, c) + g.getVal(index, c)*sVal, i, c);
	}
    }
  // this happens for spherical noise models.
  if(numCovStruct==1 && numTarget > 1)
    {
      double varSig = 0.0;
      for(int c=1; c<numTarget; c++)
	{

	  for(int i=0; i<numData; i++)
	    {
	      sVal = s.getVal(i, 0);
	      varSig = noise.getVarSigma(i, c)
		-sVal*sVal*nu.getVal(index, c);
	      noise.setVarSigma(varSig, i, c);
	      noise.setMu(noise.getMu(i, c) + g.getVal(index, c)*sVal, i, c);
	    }
	}
    }

}
int CIvm::selectPointAdd() 
{
  // returns data index of point to add.
  int index = 0;
  switch(selectionCriterion)
    {
    case RANDOM:
      index = randomPointAdd();
      break;
    case ENTROPY:
      index = entropyPointAdd();
      break;
    case RENTROPY:
      if(activeSet.size()>0)
	index = entropyPointAdd();
      else
	index = randomPointAdd();
      break;
    default:
      cerr << "Data point selection type not specified";
    }
  return index;
}
int CIvm::entropyPointAdd()
{
  // choose point from inactive set to add via entropy selection.
  vector<double> delta;
  delta.reserve(inactiveSet.size());
  for(int i=0; i<inactiveSet.size(); i++)
    delta.push_back(entropyChangeAdd(inactiveSet[i]));
  vector<double>::iterator maxVal = 
    max_element(delta.begin(), delta.end());
  changeEntropy(*maxVal);  // store global entropy change.
  return inactiveSet[maxVal - delta.begin()];
}

int CIvm::randomPointAdd() 
{
  // choose point from inactive set to add randomly.
  int index = rand();
  index = (index*inactiveSet.size())/RAND_MAX;
  index = inactiveSet[index];
  changeEntropy(entropyChangeAdd(index));
  return index;
}

double CIvm::entropyChangeAdd(const int index) const
{
  // compute the entropy change associated with point addition.
  // make sure that index is in the inactive set.
  assert(find(inactiveSet.begin(), inactiveSet.end(), index)!=inactiveSet.end());
  double entChange=0.0;
  if(noise.isSpherical())
    {
      entChange = -.5*log2(1-noise.getVarSigma(index, 0)
			   *nu.getVal(index, 0)+1e-300)*numTarget;
    }
  else
    {
      for(int j=0; j<numTarget; j++)
	entChange += -.5*log2(1-noise.getVarSigma(index, j)
			      *nu.getVal(index, j)+1e-300);
    }
  return entChange;
}
int CIvm::selectPointRemove()
{
  // returns data index of point to remove.
  int index = 0;
  switch(selectionCriterion)
    {
    case RANDOM:
      index = randomPointRemove();
      break;
    case ENTROPY:
    case RENTROPY:
      index = entropyPointRemove();
      break;
    default:
      cerr << "Data point selection type not specified";
    }
  return index;
}
int CIvm::entropyPointRemove() 
{
  vector<double> delta;
  delta.reserve(activeSet.size());
  for(int i=0; i<activeSet.size(); i++)
    delta.push_back(entropyChangeRemove(activeSet[i]));
  vector<double>::iterator maxVal = 
    max_element(delta.begin(), delta.end());
  changeEntropy(*maxVal);
  return inactiveSet[maxVal - delta.begin()];
}

int CIvm::randomPointRemove() 
{
  int index = rand();
  index = (index*activeSet.size())/RAND_MAX;
  index = activeSet[index];
  changeEntropy(entropyChangeRemove(index));
  return index;
}

double CIvm::entropyChangeRemove(const int index) const
{
  // compute entropy change associated with point removal.
  // make sure that index is in the active set.
  assert(find(activeSet.begin(), activeSet.end(), index)!=activeSet.end());
  double entChange = 0.0;
  if(noise.isSpherical())
    {
      entChange = -.5*log2(1-noise.getVarSigma(index, 0)
			   *beta.getVal(activeSet[index], 0)+1e-300)*numTarget;
    }
  else
    {
      for(int j=0; j<numTarget; j++)
	entChange += -.5*log2(1-noise.getVarSigma(index, j)
			      *beta.getVal(activeSet[index], j)+1e-300);
    }
  return entChange;
}
void CIvm::updateNuG()
{
  for(int i=0; i<numData; i++)
    noise.getNuG(g, nu, i);
}
void CIvm::updateK() const
{
  double kVal=0.0;
  for(int i=0; i<activeSet.size(); i++)
    {
      K.setVal(kern.diagComputeElement(X, activeSet[i]), i, i);
      for(int j=0; j<i; j++)
	{
	  kVal=kern.computeElement(X, activeSet[i], X, activeSet[j]);
	  K.setVal(kVal, i, j);
	  K.setVal(kVal, j, i);
	}
    }
  K.setSymmetric(true);
}
void CIvm::updateInvK(const int dim) const
{
  invK.deepCopy(K);
  for(int i=0; i<activeSetSize; i++)
    invK.setVal(invK.getVal(i, i) + 1/beta.getVal(i, dim), i, i);
  invK.setSymmetric(true);
  CMatrix U(chol(invK));
  logDetK = invK.logDet(U); 
  invK.pdinv(U);
}
  
double CIvm::approxLogLikelihood() const
{
  double L=0.0;
  updateK();
  CMatrix invKm(invK.getRows(), 1);
  if(noise.isSpherical())
    {
      updateInvK(0);
    }
  for(int j=0; j<m.getCols(); j++)
    {
      if(!noise.isSpherical())
	updateInvK(j);
      invK.setSymmetric(true);
      invKm.symvColCol(0, invK, m, j, 1.0, 0.0, "u");
      L -= .5*(logDetK + invKm.dotColCol(0, m, j));
    }
  L+=kern.priorLogProb();
  return L;
}  
void CIvm::approxLogLikelihoodGradient(CMatrix& g) const
{
  assert(g.getRows()==1);
  assert(g.getCols()==getOptNumParams());
  g.zeros();
  CMatrix tempG(1, getOptNumParams());
  updateK();
  if(noise.isSpherical())
    {
      updateInvK(0);
    }
  for(int j=0; j<m.getCols(); j++)
    {
      if(!noise.isSpherical())
	{
	  updateInvK(j);
	}
      updateCovGradient(j);
      kern.getGradTransParams(tempG, activeX, covGrad);
      g+=tempG;
      
    }
}

mxArray* CIvm::toMxArray() const
{
  int dims[1];
  dims[0]=1;
  const char* fieldNames[]={"I", "J", "m", "beta"};
  mxArray* matlabArray = mxCreateStructArray(1, dims, 4, fieldNames);

  // The I and J fields.
  vector<int> activeMatlab = activeSet;
  for(int i=0; i<activeMatlab.size(); i++)
    activeMatlab[i]++;
  mxSetField(matlabArray, 0, "I", convertMxArray(activeMatlab));
  vector<int> inactiveMatlab = inactiveSet;
  for(int i=0; i<inactiveMatlab.size(); i++)
    inactiveMatlab[i]++;
  mxSetField(matlabArray, 0, "J", convertMxArray(inactiveMatlab));
  
  // Other matrix fields.
  CMatrix tempM(numData, m.getCols());
  CMatrix tempB(numData, beta.getCols());
  for(int i=0; i<activeSet.size(); i++)
    {
      tempM.copyRowRow(activeSet[i], m, i);
      tempB.copyRowRow(activeSet[i], beta, i);
    }
  
  mxSetField(matlabArray, 0, "m", tempM.toMxArray());
  mxSetField(matlabArray, 0, "beta", tempB.toMxArray());
  return matlabArray;
}
void CIvm::fromMxArray(const mxArray* matlabArray)
{
  activeSet = mxArrayExtractVectorIntField(matlabArray, "I");
  for(int i=0; i<activeSet.size(); i++)
    activeSet[i]--;
  inactiveSet = mxArrayExtractVectorIntField(matlabArray, "J");
  for(int i=0; i<inactiveSet.size(); i++)
    inactiveSet[i]--;
  activeSetSize = activeSet.size();
  CMatrix tempM;
  tempM.fromMxArray(mxArrayExtractMxArrayField(matlabArray, "m"));
  CMatrix tempB;
  tempB.fromMxArray(mxArrayExtractMxArrayField(matlabArray, "beta"));
  m.resize(activeSetSize, tempM.getCols());
  beta.resize(activeSetSize, tempB.getCols());
  for(int i=0; i<activeSet.size(); i++)
    {
      m.copyRowRow(i, tempM, activeSet[i]);
      beta.copyRowRow(i, tempB, activeSet[i]);
    }
  assert(activeSetSize<numData);
  assert(m.getCols()==numTarget);
  assert(beta.dimensionsMatch(m));
  // TODO check that I and J cover 1:numData
}
void CIvm::optimise(const int maxIters, const int kernIters, const int noiseIters)
{
  if(kernIters>0 || noiseIters>0)
    {
      for(int iters=0; iters<maxIters; iters++)
	{
	  if(getVerbosity()>1)
	    cout << "IVM External Iteration: " << iters+1 << endl;
	  if(kernIters>0)
	    {
	      init();
	      selectPoints();
	      if(getVerbosity()>3 && getOptNumParams()<10)
		checkGradients();
	      if(getVerbosity()>1)
		cout << "Optimising kernel parameters ..." <<endl;
	      scgOptimise(kernIters);
	      if(getVerbosity()>1)
		cout << "... done. " << endl;
	      if(getVerbosity()>2)
		kern.display(cout);
	    }
	  if(noiseIters>0)
	    {
	      init();
	      selectPoints();
	      if(getVerbosity()>3 && noise.getOptNumParams()<10)
		noise.checkGradients();
	      if(getVerbosity()>1)
		cout << "Optimising noise parameters ..." <<endl;
	      noise.scgOptimise(noiseIters);
	      if(getVerbosity()>1)
		cout << "... done." <<endl;
	      if(getVerbosity()>2)
		noise.display(cout);
	    }
	}
    }
  init();
  selectPoints();
  if(getVerbosity()>0)
    display(cout);
}
bool CIvm::equals(const CIvm& model, const double tol) const
{
  if(!noise.equals(model.noise, tol))
    return false;
  if(!kern.equals(model.kern, tol))
    return false;
  if(!m.equals(model.m, tol))
    return false;
  if(!beta.equals(model.beta, tol))
    return false;
  if(activeSet!=model.activeSet)
    return false;
  if(inactiveSet!=model.inactiveSet)
    return false;
  return true;
}
void CIvm::display(ostream& os) const 
{
  cout << "IVM Model: " << endl;
  cout << "Active Set Size: " << activeSetSize << endl;
  cout << "Kernel Type: " << endl;
  kern.display(os);
  cout << "Noise Type: " << endl;
  noise.display(os);
}

void CIvm::updateCovGradient(const int index) const
{
  CMatrix invKm(invK.getRows(), 1);
  invK.setSymmetric(true);
  invKm.symvColCol(0, invK, m, index, 1.0, 0.0, "u");
  covGrad.deepCopy(invK);
  covGrad.syr(invKm, -1.0, "u");
  covGrad.scale(-0.5);
}


