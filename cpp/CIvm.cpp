#include "CIvm.h"

CIvm::CIvm(const CMatrix& inData, const CMatrix& targetData, 
	   CKern& kernel, CNoise& noiseModel, const int selectCrit,
	   const int dVal) 
  : X(inData), y(targetData), kern(kernel), 
    noise(noiseModel), selectionCriterion(selectCrit), numTarget(y.getCols()), numData(y.getRows()), lastEntropyChange(0.0), cumEntropy(0.0), activeSetSize(dVal)
{
  terminate = false;
  epUpdate = false;
  // need to set up Sigma for handing spherical and non-spherical.
  // set numCovStruct
  if(noise.isSpherical())
    numCovStruct = 1;
  else
    numCovStruct = numTarget;
  init();
}
void CIvm::init()
{
  // set Kstore to zeros numData, activeSetSize.
  Kstore.resize(numData, activeSetSize);
  Kstore.setVals(0.0);
  // set m and beta to zeros(size of y) -- could do with sparse representation.
  m.resize(activeSetSize, y.getCols());
  beta.resize(activeSetSize, y.getCols());
  // set nu to zeros(size of y)
  nu.resize(y.getRows(), y.getCols());
  nu.setVals(0.0);
  // set g to zeros(size of y)
  g.resize(y.getRows(), y.getRows());
  g.setVals(0.0);
  // set mu and varSigma to zeros(size of y)
  mu.resize(y.getRows(), y.getCols());
  mu.setVals(0.0);
  varSigma.resize(y.getRows(), y.getCols());
  double dk=0.0;
  for(int i=0; i<varSigma.getRows(); i++)
    {
      dk = kern.diagComputeElement(X, i);
      for(int j=0; j<varSigma.getCols(); j++)
	{
	  varSigma.setVals(dk, i, j);
	}
    }

  // set up K invK and covGrad
  K.resize(activeSetSize, activeSetSize);
  invK.resize(activeSetSize, activeSetSize);
  covGrad.resize(activeSetSize, activeSetSize);
  activeX.resize(activeSetSize, X.getCols());
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
      M[c].setVals(0.0);
      L[c].resize(activeSetSize, activeSetSize);
      L[c].setVals(0.0);
      Linv[c].resize(activeSetSize, activeSetSize);
      Linv[c].setVals(0.0);
    }

  // set up covGrad and invK storage matrices
  invK.resize(activeSetSize, activeSetSize);
  invK.zeros();
  invK.setSymmetric(true);
  covGrad.resize(activeSetSize, activeSetSize);
  covGrad.zeros();
  covGrad.setSymmetric(true);

  // fill the inactive set.
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
  for(int k=0; k<activeSetSize; k++)
    {
      index = selectPointAdd();
      addPoint(index);
      
    }
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
  noise.updateSites(m, beta, actIndex, g, nu, mu, varSigma, y, index);
  for(int j=0; j<beta.getCols(); j++)
    {
      double betVal = beta.getVals(actIndex, j);
      if(betVal<0)
	{
	  if(noise.isLogConcave())
	    {
	      cerr << "Error beta less than zero for log concave model.";
	    }
	  else
	    {
	      beta.setVals(1e-6, actIndex, j);
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
      Kstore.setVals(kern.computeElement(X, i, X, index), i, activePoint);      
    }
  // add white noise term to relevant index.
  double val = Kstore.getVals(index, activePoint);
  Kstore.setVals(val+kern.getWhite(), index, activePoint);
  double lValInv = 0.0;
  double vs = 0.0;
  double ms = 0.0;
  double sVal = 0.0;
  for(int c=0; c<numCovStruct; c++)
    {
      assert(nu.getVals(index, c)>=0);
      lValInv = sqrt(nu.getVals(index, c));
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
      L[c].setVals(1/lValInv, activePoint, activePoint);
      a.trans(); // turn a into a column vector.
      // update the varSigma and mu systems.
      for(int i=0; i<numData; i++)
	{
	  sVal = s.getVals(i, 0);
	  varSigma.setVals(varSigma.getVals(i, c)
			   -sVal*sVal*nu.getVals(index, c), i, c);
	  assert(varSigma.getVals(i, c)>=0);
	  mu.setVals(mu.getVals(i, c) + g.getVals(index, c)*sVal, i, c);
	}
    }
  // this happens for spherical noise models.
  if(numCovStruct==1 && numTarget > 1)
    {
      for(int c=1; c<numTarget; c++)
	{
	  for(int i=0; i<numData; i++)
	    {
	      sVal = s.getVals(i, 0);
	      varSigma.setVals(varSigma.getVals(i, c)
			       -sVal*sVal*nu.getVals(index, c), i, c);
	      mu.setVals(mu.getVals(i, c) + g.getVals(index, c)*sVal, i, c);
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
      entChange = -.5*log2(1-varSigma.getVals(index, 0)
			   *nu.getVals(index, 0)+1e-300)*numTarget;
    }
  else
    {
      for(int j=0; j<numTarget; j++)
	entChange += -.5*log2(1-varSigma.getVals(index, j)
			      *nu.getVals(index, j)+1e-300);
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
      entChange = -.5*log2(1-varSigma.getVals(index, 0)
			   *beta.getVals(activeSet[index], 0)+1e-300)*numTarget;
    }
  else
    {
      for(int j=0; j<numTarget; j++)
	entChange += -.5*log2(1-varSigma.getVals(index, j)
			      *beta.getVals(activeSet[index], j)+1e-300);
    }
  return entChange;
}
void CIvm::updateNuG()
{
  for(int i=0; i<numData; i++)
    noise.getNuG(g, nu, mu, varSigma, y, i);
}
void CIvm::updateK()
{
  double kVal=0.0;
  for(int i=0; i<activeSet.size(); i++)
    {
      K.setVals(kern.diagComputeElement(X, activeSet[i]), i, i);
      for(int j=0; j<i; j++)
	{
	  kVal=kern.computeElement(X, activeSet[i], X, activeSet[j]);
	  K.setVals(kVal, i, j);
	  K.setVals(kVal, j, i);
	}
    }
  K.setSymmetric(true);
}
void CIvm::updateInvK(const int dim)
{
  invK.deepCopy(K);
  for(int i=0; i<activeSetSize; i++)
    invK.setVals(invK.getVals(i, i) + 1/beta.getVals(i, dim), i, i);
  invK.setSymmetric(true);
  invK.writeMatlabFile("invK.mat", "K");
  //TODO this should be done in tandem (log det and invK).
  CMatrix U(chol(K));
  logDetK = 0.0;
  for(int i=0; i<U.getRows(); i++)
    logDetK+=log(U.getVals(i, i));
  logDetK *= 2;
  invK.pdinv();
}
  
double CIvm::approxLogLikelihood()
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
      L -= .5*(logDetK + invKm.dotCol(0, m, j));
    }
  return L;
}  
void CIvm::approxLogLikelihoodGradient(CMatrix& g)
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
 
void CIvm::updateCovGradient(const int index)
{
  CMatrix invKm(invK.getRows(), 1);
  invK.setSymmetric(true);
  invKm.symvColCol(0, invK, m, index, 1.0, 0.0, "u");
  covGrad.deepCopy(invK);
  covGrad.syr(invKm, -1.0, "u");
  covGrad.scale(-0.5);
}


