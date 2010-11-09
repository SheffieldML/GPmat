#include "CGp.h"
CGp::CGp() 
  : CMapModel(), CProbabilisticOptimisable()
{
  _init();
}

CGp::CGp(unsigned int q, unsigned int d, 
	 CMatrix* pXin, CMatrix* pyin, 
	 CKern* pkernel, CNoise* pnois, 
	 int approxType, 
	 unsigned int actSetSize, int verbos)
  : CMapModel(), CProbabilisticOptimisable(), pkern(pkernel), pnoise(pnois), pX(pXin), py(pyin)
{
  DIMENSIONMATCH(pXin->getCols()==q);
  DIMENSIONMATCH(pyin->getCols()==d);
  DIMENSIONMATCH(pyin->getRows()==pXin->getRows());

  pnoise->setTarget(pyin);

  _init();
  numActive = actSetSize;
  setApproximationType(approxType);
  if(isSparseApproximation() && numActive<1)
    throw ndlexceptions::Error("Error sparse approximations need active sets of non-zero size.");

  setVerbosity(verbos);
  initKernStoreage();

  setOutputDim(pnoise->getOutputDim());
  setInputDim(pXin->getCols());
  setNumData(pnoise->getNumData());
  initStoreage();
  // set scale and bias.
  scale.setVals(1.0);
  bias.zeros();
  initVals();
  
  switch(approxType) 
  {
  case FTC:
  case DTC:
  case DTCVAR:
  case FITC:
  case PITC:
    // positive constraint on beta
    betaTransform = CTransform::defaultPositive();
    break;
  default:
    throw ndlexceptions::Error("Unrecognised sparse approximation type.");
  }
}
CGp::CGp(CKern* pkernel, CNoise* pnois,
	 CMatrix* pXin, int approxType, 
	 unsigned int actSetSize, int verbos)
  :
  CMapModel(), CProbabilisticOptimisable(),
  pkern(pkernel), pnoise(pnois), pX(pXin)
{
  DIMENSIONMATCH(pXin->getRows()==pnoise->getNumData());
  _init();
  numActive = actSetSize;
  setApproximationType(approxType);
  if(isSparseApproximation() && numActive<1)
    throw ndlexceptions::Error("Error sparse approximations need active sets of non-zero size.");
  setVerbosity(verbos);
  initKernStoreage();

  setOutputDim(pnoise->getOutputDim());
  setInputDim(pXin->getCols());
  setNumData(pnoise->getNumData());
  py = pnoise->py;
  initStoreage();
  // set scale and bias.
  scale.setVals(1.0);
  bias.zeros();
  initVals();
  
  switch(approxType) 
  {
  case FTC:
  case DTC:
  case DTCVAR:
  case FITC:
  case PITC:
    // positive constraint on beta
    betaTransform = CTransform::defaultPositive();
    break;
  default:
    throw ndlexceptions::Error("Unrecognised sparse approximation type.");
  }
}
void CGp::_init()
{
  setType("gp");
  setName("Gaussian process");
  setOutputScaleLearnt(false);
  setOutputBiasLearnt(false);
  jitter = 1e-6;
  inducingFixed=false;
  spherical=true;
  optimiseX = false;
  backConstrained = false;
  KupToDate = false;
  sparseApproximation = false;
  setApproximationType(FTC);
  numActive = 0;
}

void CGp::initKernStoreage() 
{
  g_param.resize(1, pkern->getNumParams());
}

void CGp::initOptimiseXStoreage()
{
  if(isOptimiseX())
  {
    if(!isBackConstrained())
    {
      if(!isSparseApproximation())
      {
	dgKX.resize(getNumData(), getInputDim()); /// TODO: this also may be too big.
	gKX.resize(getNumData(), getInputDim());
      }
      else
      {
	gKX_uf2.resize(numActive, getInputDim());
      }
      gXorW.resize(getNumData(), getInputDim());
    }
    else
    {
      gXorW.resize(1, pbackModel->getOptNumParams());
    }
  }
}
void CGp::initSparseStoreage()
{
  if (isSparseApproximation()) 
  {
    X_u.resize(numActive, getInputDim());
    if (!isInducingFixed()) 
    {
      gX_u.resize(numActive, getInputDim());
    }
    E.resize(numActive, getOutputDim());
    EET.resize(numActive, numActive);
    AinvEET.resize(numActive, numActive);
    AinvEETAinv.resize(numActive, numActive);
    Ainv.resize(numActive, numActive);
    AinvK_uf.resize(numActive, getNumData());
    EMT.resize(numActive, getNumData());
    AinvEMT.resize(numActive, getNumData());
    A.resize(numActive, numActive);
    K_uu.resize(numActive, numActive);
    invK_uu.resize(numActive, numActive);
    K_uf.resize(numActive, getNumData());
    LcholK.resize(numActive, numActive);
    gK_uf.resize(numActive, getNumData());
    gK_uu.resize(numActive, numActive);
    gKX.resize(numActive, getInputDim());
    dgKX.resize(numActive, getInputDim());
    gKX_uf.resize(getNumData(), getInputDim());
    gBeta.resize(1, 1);
    Alpha.resize(numActive, getOutputDim());
  }
  else 
  {
    K.resize(getNumData(), getNumData());
    invK.resize(getNumData(), getNumData());
    LcholK.resize(getNumData(), getNumData());    
    covGrad.resize(getNumData(), getNumData());  
    Alpha.resize(getNumData(), getOutputDim());
  }
}
void CGp::initStoreage()
{
  //py->resize(getNumData(), getOutputDim());
  m.resize(getNumData(), getOutputDim());
  setMupToDate(false);

  // TODO cater for all different types of beta.
  if(isSpherical()) 
  {
    // TODO : This should only be one value! But the noise initialisation requires large matrix
    beta.resize(getNumData(), getOutputDim());
  }
  else 
  {
    beta.resize(getNumData(), getOutputDim());
  }
  
  scale.resize(1, getOutputDim());
  bias.resize(1, getOutputDim());
  
  // gradient matrices
  if (isOutputScaleLearnt()) 
  {
    g_scaleBias.resize(1, getOutputDim());
  }
  nu.resize(getNumData(), getOutputDim());
  g.resize(getNumData(), getOutputDim());

  initOptimiseXStoreage();
  initSparseStoreage();
   
  // approximation specific allocations
  switch(getApproximationType()) 
  {
  case FTC:
    break;
  case DTC:
    break;
  case DTCVAR:
    gLambda.resize(getNumData(), 1);
    diagK.resize(getNumData(), 1);
    invK_uuK_uf.resize(numActive, getNumData());
    V.resize(numActive, getNumData());  
    break;
  case FITC:
    // Temporary variables for FITC approximation.
    gLambda.resize(getNumData(), 1);
    diagK.resize(getNumData(), 1);
    scaledM.resize(getNumData(), getOutputDim());
    V.resize(numActive, getNumData());  
    Am.resize(numActive, numActive);  
    Lm.resize(numActive, numActive);  
    invLmV.resize(numActive, getNumData());  
    bet.resize(numActive, getOutputDim());
    AinvE.resize(numActive, getOutputDim());
    diagMMT.resize(getNumData(), 1);
    diagQ.resize(getNumData(), 1);
    diagK_ufdAinvplusAinvEETAinvK_fu.resize(getNumData(), 1); 
    invK_uuK_uf.resize(numActive, getNumData());
    invK_uuK_ufDinv.resize(numActive, getNumData());
    invK_uuK_ufDinvQ.resize(numActive, getNumData());
    K_ufdotTimesAinvEMT.resize(numActive, getNumData());
    diagK_ufAinvEMT.resize(getNumData(), 1);
    break;
  case PITC:
    break;
  }    
   
}

void CGp::updateM() const 
{
  if(!isMupToDate()) 
  {
    for(unsigned int i=0; i<getOutputDim(); i++) 
    {
      m.copyColCol(i, *py, i);
      m.addCol(i, -bias.getVal(i));
      m.scaleCol(i, 1/scale.getVal(i));
    }
    setMupToDate(true);
  }
}
void CGp::initVals() 
{
  setSpherical(true); // not implemented non-spherical stuff yet.
  maxTries = 10; // number of cholesky decompositions to be attempted (adding jitter each time)
  // TODO: this isn't being properly used ...
  for(size_t i=0; i<getNumData(); i++) 
  {
    pnoise->updateSites(m, beta, i, nu, g, i);
  }
  setBetaVals(1e3);
  updateM();
  
  if(isSparseApproximation()) 
  {
    // Randomise selection of active points.
    if(numActive>getNumData())
      throw ndlexceptions::CommandLineError("Number of active points has to be less than number of data.");
    vector<unsigned long> ind = ndlutil::randpermTrunc(getNumData(), numActive);
    sort(ind.begin(), ind.end());
    for(unsigned int i=0; i<ind.size(); i++) 
    {
      X_u.copyRowRow(i, *pX, ind[i]);	
    }
  }
}


void CGp::updateX()
{
  setKupToDate(false);
}
string CGp::getNoiseType() const
{
  return pnoise->getType();
}

unsigned int CGp::getOptNumParams() const
{
  int tot = pkern->getNumParams();
  
  // for GP-LVM
  if(isOptimiseX())
  {
    if(!isBackConstrained())
    {
      tot += getInputDim()*getNumData();
    }
    else
    {
      throw ndlexceptions::NotImplementedError("Back constraints not yet implemented.");
    }
  }
  // for FITC, DTC, DTCVAR and PITC
  if(isSparseApproximation()) 
  {
    if(!isInducingFixed()) 
    {
      tot += numActive*getInputDim();
    }
    tot++; // for beta value.
  }
  // for learning scales.
  if(isOutputScaleLearnt()) 
  {
    tot+=getOutputDim();
  }
  return tot;
}

void CGp::getOptParams(CMatrix& param) const
{
  int counter = 0;
  // for GP-LVM
  if(isOptimiseX())
  {
    if(!isBackConstrained())
    {
      for(unsigned int j=0; j<getInputDim(); j++)
      {
	for(unsigned int i=0; i<getNumData(); i++)
	{
	  param.setVal(pX->getVal(i, j), counter);
	  counter++;
	}
      }
    }
    else
    {
      throw ndlexceptions::NotImplementedError("Back constraints not yet implemented.");
    }
  }
  // for FITC, DTC, DTCVAR and PITC
  if(isSparseApproximation()) 
  {
    if(!isInducingFixed()) 
    {
      for(unsigned int j=0; j<getInputDim(); j++) 
      { 
	for(unsigned int i=0; i<numActive; i++) 
	{
	  param.setVal(X_u.getVal(i, j), counter);
	  counter++;
	}
      }
    }
  }
  for(unsigned int i=0; i<pkern->getNumParams(); i++)
  {
    param.setVal(pkern->getTransParam(i), counter);
    counter++;
  }
  if(isOutputScaleLearnt()) 
  {
    for(unsigned int j=0; j<getOutputDim(); j++) 
    {
      param.setVal(getScaleVal(j), counter);
      counter++;
    }
  }
  if(isSparseApproximation()) 
  {
    param.setVal(betaTransform->xtoa(beta.getVal(0, 0)), counter);
    counter++;
  }
}

void CGp::setOptParams(const CMatrix& param)
{
  setKupToDate(false);
  int counter=0;
  if(isOptimiseX())
  {
    if(!isBackConstrained())
    {
      for(unsigned int j=0; j<getInputDim(); j++)
      {
	for(unsigned int i=0; i<getNumData(); i++)
	{
	  pX->setVal(param.getVal(counter), i, j);
	  counter++;
	}
      }
    }
    else
    {
      throw ndlexceptions::NotImplementedError("Back constraints not yet implemented.");
    }
  }
  if(isSparseApproximation()) 
  {
    if(!isInducingFixed()) 
    {
      for(unsigned int j=0; j<getInputDim(); j++) 
      {
        for(unsigned int i=0; i<numActive; i++) 
        {
          X_u.setVal(param.getVal(counter), i, j);
          counter++;
        }
      }
    }
  }
  for(unsigned int i=0; i<pkern->getNumParams(); i++) 
  {	  
    pkern->setTransParam(param.getVal(counter), i);
    counter++;
  }
  if(isOutputScaleLearnt())
  {
    for(unsigned int j=0; j<getOutputDim(); j++) 
    {
      scale.setVal(param.getVal(counter), j);
      counter++;
    }      
    updateM();
  }
  if(isSparseApproximation()) 
  {
    beta.setVal(betaTransform->atox(param.getVal(counter)), 0);
    counter++;
  }
  
}

void CGp::out(CMatrix& yPred, const CMatrix& Xin) const
{
  DIMENSIONMATCH(yPred.getRows()==Xin.getRows());
  CMatrix muTest(yPred.getRows(), yPred.getCols());
  CMatrix varSigmaTest(yPred.getRows(), yPred.getCols());
  posteriorMeanVar(muTest, varSigmaTest, Xin);
  pnoise->out(yPred, muTest, varSigmaTest);
}

void CGp::out(CMatrix& yPred, CMatrix& probPred, const CMatrix& Xin) const
{
  CMatrix muTest(yPred.getRows(), yPred.getCols());
  CMatrix varSigmaTest(yPred.getRows(), yPred.getCols());
  posteriorMeanVar(muTest, varSigmaTest, Xin);
  pnoise->out(yPred, probPred, muTest, varSigmaTest);
}
double CGp::outGradParams(CMatrix& g, const CMatrix &Xin, unsigned int pointNo, unsigned int outputNo) const
{
  throw ndlexceptions::NotImplementedError("outGradParams not yet implemented for CGp.");
}
double CGp::outGradX(CMatrix& g, const CMatrix &Xin, unsigned int pointNo, unsigned int outputNo) const
{
  throw ndlexceptions::NotImplementedError("outGradX not yet implemented for CGp.");
}
void CGp::updateAlpha() const
{ 
  if(!isAlphaUpToDate())
  { 
    updateM();
    updateK();
    updateAD();  
    switch(getApproximationType())
    {
    case FTC:
      if(isSpherical())
      {
	Alpha.deepCopy(m);
	Alpha.trsm(LcholK, 1.0, "l", "l", "n", "n");
	Alpha.trsm(LcholK, 1.0, "l", "l", "t", "n");
      }
      else
      {
	throw ndlexceptions::Error("Non-spherical implementations not yet in place for FTC");
      }
      break;
    case DTC:
    case DTCVAR:
      if(isSpherical())
      {
	Alpha.gemm(K_uf, m, 1.0, 0.0, "n", "n");
	Alpha.trsm(LcholA, 1.0, "l", "l", "n", "n");
	Alpha.trsm(LcholA, 1.0, "l", "l", "t", "n");
      }
      else
      {
	throw ndlexceptions::Error("Non-spherical implementations not yet in place for DTC");
	
      }
      break;
    case FITC:
      if(isSpherical())
      {
	//////////////
	scaledM.deepCopy(m);
	for(unsigned int i=0; i<getNumData(); i++)
	  scaledM.scaleRow(i, 1/diagD.getVal(i));
	/////////////
	setADupToDate(false);
	Alpha.gemm(K_uf, scaledM, 1.0, 0.0, "n", "n");
	Alpha.trsm(LcholA, 1.0, "l", "l", "n", "n");
	Alpha.trsm(LcholA, 1.0, "l", "l", "t", "n");
      }
      else
      {
        throw ndlexceptions::Error("Non-spherical implementations not yet in place for FITC");
      }
      break;
    case PITC:
      if (isSpherical()) 
      {
      }
      else 
      {
        throw ndlexceptions::Error("Non-spherical implementations not yet in place for PITC");
      }
      break;
    }
    setAlphaUpToDate(true);
  }
}
void CGp::_testComputeKx(CMatrix &kX, const CMatrix& Xin) const
{
  if(isSparseApproximation())
  {
    DIMENSIONMATCH(kX.getRows()==X_u.getRows() && kX.getCols()==Xin.getRows());  
    pkern->compute(kX, X_u, Xin);
  }
  else
  {
    DIMENSIONMATCH(kX.getRows()==pX->getRows() && kX.getCols()==Xin.getRows());  
    pkern->compute(kX, *pX, Xin);
  }
}
void CGp::_posteriorMean(CMatrix& mu, const CMatrix& kX) const
{
  updateAlpha();
  DIMENSIONMATCH(mu.getRows()==kX.getCols());
  DIMENSIONMATCH(mu.getCols()==getOutputDim());
  for(unsigned int i=0; i<kX.getCols(); i++)
  {
    for(unsigned int j=0; j<getOutputDim(); j++)
    {
      mu.setVal(Alpha.dotColCol(j, kX, i), i, j);
    }
  }
  // apply output scale and bias.
  for(unsigned int j=0; j<getOutputDim(); j++)
  {
    double scaleVal = scale.getVal(j);
    if(scaleVal!=1.0)
    {
      mu.scaleCol(j, scaleVal);
    }
    double biasVal = bias.getVal(j);
    if(biasVal!=0.0)
    {
      mu.addCol(j, biasVal);
    }
  }
}
void CGp::_posteriorVar(CMatrix& varSigma, CMatrix& kX, const CMatrix& Xin) const
{
  updateK();
  updateAD();
  DIMENSIONMATCH(varSigma.getRows()==kX.getCols());
  DIMENSIONMATCH(varSigma.getCols()==getOutputDim());
  // WARNING: destroys kX through in place operations.
  if(isSparseApproximation())
  {
    CMatrix store(kX.getRows(), kX.getCols());
    EET.deepCopy(invK_uu);
    EET.axpy(Ainv, -1.0/getBetaVal());
    EET.setSymmetric(true);
    store.symm(EET, kX, 1.0, 0.0, "l", "l"); 
    for(unsigned int i=0; i<kX.getCols(); i++)
    {
      double vsVal = pkern->diagComputeElement(Xin, i) - kX.dotColCol(i, store, i);
      CHECKZEROORPOSITIVE(vsVal>=0);
      vsVal += 1.0/getBetaVal();
      for(unsigned int j=0; j<getOutputDim(); j++)
      {
	varSigma.setVal(vsVal, i, j);
      }
    }
  }
  else 
  {
    // solve LcholK * tmp = kX,  kX := tmp
    kX.trsm(LcholK, 1.0, "L", "L", "N", "N"); // now it is LcholKinv*K
    for(unsigned int i=0; i<kX.getCols(); i++)
    {
      double vsVal = pkern->diagComputeElement(Xin, i) - kX.norm2Col(i);
      CHECKZEROORPOSITIVE(vsVal>=0);
      for(unsigned int j=0; j<getOutputDim(); j++)
      {
        varSigma.setVal(vsVal, i, j);
      }
    }
  }
  
  
  // apply output scale and bias.
  for(unsigned int j=0; j<getOutputDim(); j++)
  {
    double scaleVal = scale.getVal(j);
    if(scaleVal!=1.0)
    {
      varSigma.scaleCol(j, scaleVal*scaleVal);
    }
  }
}
void CGp::posteriorMean(CMatrix& mu, const CMatrix& Xin) const
{
  updateAlpha();
  int alRows = 0;
  if(isSparseApproximation())
  {
    alRows = numActive;
  }
  else
  { 
    alRows = getNumData();
  }
  CMatrix kX(alRows, Xin.getRows());
  _testComputeKx(kX, Xin);
  _posteriorMean(mu, kX);
}
void CGp::posteriorMeanVar(CMatrix& mu, CMatrix& varSigma, const CMatrix& Xin) const
{
  DIMENSIONMATCH(mu.getCols()==getOutputDim());
  DIMENSIONMATCH(varSigma.getCols()==getOutputDim());
  DIMENSIONMATCH(mu.getRows()==Xin.getRows());
  DIMENSIONMATCH(varSigma.getRows()==Xin.getRows());
  
  int alRows = 0;
  if(isSparseApproximation())
  {
    alRows = numActive;
  }
  else
  { 
    alRows = getNumData();
  }
  CMatrix kX(alRows, Xin.getRows());
  _testComputeKx(kX, Xin);
  _posteriorMean(mu, kX);
  _posteriorVar(varSigma, kX, Xin); // destroys kX through in place operations.
  
}

// Gradient routines
void CGp::updateCovGradient(unsigned int j, CMatrix& invKm) const
{
  // covGrad := 0.5*(invK * Y(:,j) * Y^t(:,j) * invK - invK)
  
  // (The 'D' factor doesn't need to be on the last term because we add 
  //  D of these covGrad expressions together in a loop in the caller.  
  //  In other words, this is covGrad for a single process only.)
  
  invKm.resize(invK.getRows(), 1);
  invKm.symvColCol(0, invK, m, j, 1.0, 0.0, "u"); // invKm = invK*m(:,index)
  covGrad.deepCopy(invK);                             // covgrad = invK
  covGrad.syr(invKm, -1.0, "u");                 // covGrad -= invKm * invKm^t
  covGrad.scale(-0.5);                           // covGrad *= -0.5
}


void CGp::updateK() const
{
  if(!isKupToDate())
  {
    _updateK();
    _updateInvK();
    setADupToDate(false);
    setKupToDate(true);
  }
}

void CGp::_updateK() const
{
  double kVal=0.0;
  switch(getApproximationType()) 
  {
  case FTC:
    // TODO: These computes should be done with pkern->compute, which 
    // could be made multi-threaded.
    for(unsigned int i=0; i<getNumData(); i++) 
    {
      K.setVal(pkern->diagComputeElement(*pX, i), i, i);
      for(unsigned int j=0; j<i; j++) 
      {
	kVal=pkern->computeElement(*pX, i, *pX, j);
	K.setVal(kVal, i, j);
	K.setVal(kVal, j, i);
      }
    }
    K.setSymmetric(true);
    break;
  case DTC:
  case DTCVAR:
  case FITC:
  case PITC:
    // TODO: These computes should be done with pkern->compute, which 
    // could be made multi-threaded.
    for(unsigned int i=0; i<numActive; i++) 
    {
      K_uu.setVal(pkern->diagComputeElement(X_u, i), i, i);
      for(unsigned int j=0; j<i; j++) 
      {
	kVal=pkern->computeElement(X_u, i, X_u, j);
	K_uu.setVal(kVal, i, j);
	K_uu.setVal(kVal, j, i);
      }
      K_uu.setSymmetric(true);
      for(unsigned int j=0; j<getNumData(); j++) {
	kVal=pkern->computeElement(X_u, i, *pX, j);
	K_uf.setVal(kVal, i, j);
      }
    }
    break;
  }
  switch(getApproximationType()) 
  {
  case FTC:
    break;
  case DTC:
    break;
  case FITC:
  case DTCVAR:
    for(unsigned int i=0; i< getNumData(); i++) 
    {
      diagK.setVal(pkern->diagComputeElement(*pX, i), i);
    }
    break;
  case PITC:
    break;
  default:
    break;
  }
  
}

void CGp::updateAD() const {
  if(!isADupToDate())
  {
    double betaVal = beta.getVal(0);
    switch(getApproximationType()) 
    {
    case FTC:
      // TODO update the inner products here.
      break;
    case DTC:
    case DTCVAR:
      if (isSpherical()) 
      {
        A.deepCopy(K_uu);
        A.gemm(K_uf, K_uf, 1.0, 1.0/betaVal, "n", "t");
        A.setSymmetric(true);
	double jit = LcholA.jitChol(A);
	if(jit>1e-2)
	  if(getVerbosity()>2)
	    cout << "Warning: jitter of " << jit << " added to A in updateAD()." << endl;

	logDetA = logDet(LcholA);
	Ainv.setSymmetric(true);
        Ainv.pdinv(LcholA);
	// Now stored as a lower cholesky.
        LcholA.trans();
	if(getApproximationType()==DTCVAR)
	{
	  V.gemm(invK_uu, K_uf, 1.0, 0.0, "n", "n");
	  V.dotMultiply(K_uf);

	  diagD.deepCopy(diagK);
	  diagD.sumCol(V, -1.0, 1.0);
	  diagD.scale(getBetaVal());
	}	
      }
      else 
      {
        throw ndlexceptions::NotImplementedError("Non-spherical implementations not yet in place for DTC or DTCVAR");
      }
      // TODO generate the inner products here
      break;
    case FITC:  
      if (isSpherical()) 
      {
// 	V.deepCopy(K_uf);
// 	V.trsm(LcholK, 1.0, "l", "l", "n", "n");
// 	V.trsm(LcholK, 1.0, "l", "l", "t", "n");
// 	V.dotMultiply(K_uf);
	V.gemm(invK_uu, K_uf, 1.0, 0.0, "n", "n");
	V.dotMultiply(K_uf);
	diagD.deepCopy(diagK);
	diagD.negate();
	diagD.sumCol(V, 1.0, 1.0);
	diagD.scale(getBetaVal());
	diagD.negate();
	diagD.add(1.0);
	
	V.deepCopy(K_uf);
	double d = 0.0;
	scaledM.deepCopy(m);
	for(unsigned int j=0; j<getNumData(); j++) 
	{
	  d = 1/diagD.getVal(j);
	  V.scaleCol(j, d);
	  scaledM.scaleRow(j, sqrt(d));
	}
	A.deepCopy(K_uu);
	A.gemm(K_uf, V, 1.0, 1.0/getBetaVal(), "n", "t");
	A.setSymmetric(true);
	// This is initially an upper Cholesky.
	double jit = LcholA.jitChol(A);
	if(jit>1e-2)
	  if(getVerbosity()>2)
	    cout << "Warning: jitter of " << jit << " added to A in updateAD()." << endl;

	logDetA = logDet(LcholA);
	Ainv.setSymmetric(true);
	Ainv.pdinv(LcholA);
	// now make it a lower Cholesky.
	LcholA.trans();
	
	V.deepCopy(K_uf);
	V.trsm(LcholK, 1.0, "l", "l", "n", "n");
	for(unsigned int j=0; j<getNumData(); j++)
	  V.scaleCol(j, 1/sqrt(diagD.getVal(j)));
	Am.diag(1/getBetaVal());
	Am.gemm(V, V, 1.0, 1.0, "n", "t");
	Am.setSymmetric(true);

	jit = Lm.jitChol(Am); // this will initially be upper triangular.
	if(jit>1e-2)
	  if(getVerbosity()>2)
	    cout << "Warning: jitter of " << jit << " added to Am in updateAD()." << endl;
	Lm.trans(); // now it is lower

	invLmV.deepCopy(V);
	invLmV.trsm(Lm, 1.0, "l", "l", "n", "n");
	bet.gemm(invLmV, scaledM, 1.0, 0.0, "n", "n");
      }
      else 
      {
	throw ndlexceptions::NotImplementedError("Non-spherical implementations not yet in place for FITC");
      }
      break;
    case PITC:
      if (isSpherical()) 
      {
	throw ndlexceptions::NotImplementedError("PITC approximation not yet implemented.");
      }
      else 
      {
	throw ndlexceptions::NotImplementedError("Non-spherical implementations not yet in place for PITC");
      }
      break;
    }
    setADupToDate(true);
  }
}
// update invK with the inverse of the kernel plus beta terms computed from the active points.
void CGp::_updateInvK(unsigned int dim) const
{
  double jit = 0.0;
  switch(getApproximationType()) {
  case FTC:
    jit = LcholK.jitChol(K); // this will initially be upper triangular.
    if(jit>1e-2)
      if(getVerbosity()>2)
	cout << "Warning: jitter of " << jit << " added to K in _updateInvK()." << endl;

    logDetK = logDet(LcholK);
    invK.setSymmetric(true);
    invK.pdinv(LcholK);
    LcholK.trans();
    break;
  case DTC:
  case DTCVAR:
  case FITC:
  case PITC:
    jit = LcholK.jitChol(K_uu); // this will initially be upper triangular.
    if(jit>1e-2)
      if(getVerbosity()>2)
	cout << "Warning: jitter of " << jit << " added to K_uu in _updateInvK()." << endl;


    logDetK_uu = logDet(LcholK);
    invK_uu.setSymmetric(true);
    invK_uu.pdinv(LcholK);
    // make it lower
    LcholK.trans();
    break;
  }
}


// compute the approximation to the log likelihood.
double CGp::logLikelihood() const
{
  updateK();
  updateAD();
  double L=0.0;
  switch(getApproximationType()) 
  {
  case FTC:
    if(isSpherical()) 
    {
      CMatrix invKm(invK.getRows(), 1);
      for(unsigned int j=0; j<getOutputDim(); j++) 
      {
	// This computes trace(invK*M*M'), column by column of M
	// invKm := invK * m(:,j)
	invKm.symvColCol(0, invK, m, j, 1.0, 0.0, "u");
	// L += invKm' * m(:,j)
	L += invKm.dotColCol(0, m, j);
	L += logDetK; 
      }
    }
    else 
    {
      throw ndlexceptions::Error("Non-spherical implementations not yet in place");
    }
    break;
  case DTC:    
  case DTCVAR:
    if(isSpherical()) 
    {
      CMatrix e(numActive, 1);
      CMatrix invAe(numActive, 1);
      L+= (double)getOutputDim()*
      (((double)numActive-(double)getNumData())*log(getBetaVal())
       -logDetK_uu + logDetA);
      for(unsigned int j=0; j<getOutputDim(); j++) 
      {
	e.gemvColCol(0, K_uf, m, j, 1.0, 0.0, "n");
	invAe.symvColCol(0, Ainv, e, 0, 1.0, 0.0, "l");
	L -= getBetaVal()*(invAe.dotColCol(0, e, 0)-m.dotColCol(j, m, j));
      }	  
      if(getApproximationType()==DTCVAR)
	L += (double)getOutputDim()*diagD.sum();
    }
    else 
    {
      throw ndlexceptions::Error("Non-spherical implementations not yet in place");
    }
    break;
    
  case FITC:
    if(isSpherical()) 
    {
      L+=((double)numActive - (double)getNumData())*log(getBetaVal()) + (double)getNumData()*ndlutil::LOGTWOPI;
      for(unsigned int i=0; i<getNumData(); i++) 
      {
	L += log(diagD.getVal(i));
      }
      double temp = 0.0;
      for(unsigned int i=0; i<numActive; i++) 
      {
	temp += log(Lm.getVal(i, i));
      }
      L += temp*2.0;
      L *= (double)getOutputDim();
      for(unsigned int j=0; j<getOutputDim(); j++) 
      {
	L+= getBetaVal()*(scaledM.dotColCol(j, scaledM, j)-bet.dotColCol(j, bet, j));
	
      }
    }
    else 
    {
      throw ndlexceptions::NotImplementedError("Non-spherical implementations not yet in place");
    }
    break;

  case PITC:
    if(isSpherical()) 
    {
      throw ndlexceptions::NotImplementedError("PITC not yet implemented.");
    }
    else 
    {
      throw ndlexceptions::NotImplementedError("Non-spherical implementations not yet in place");
    }
    break;
    
  }
  if(isOutputScaleLearnt()) 
  {
    // scales lead to terms of the form log w_j to be added.
    for(unsigned int j=0; j<getOutputDim(); j++) 
    {
      L+=2*log(fabs(scale.getVal(j)));
    }
  }
  L*=-0.5;
  L+=pkern->priorLogProb();
  L-=(double)getOutputDim()*(double)getNumData()*ndlutil::HALFLOGTWOPI;
  return L;
}
// compute the gradients wrt parameters and latent variables.
double CGp::logLikelihoodGradient(CMatrix& g) const 
{
  updateG();
  g.zeros();
  int counter = 0;
  if(isOptimiseX())
  {
    if(!isBackConstrained())
    {
      for(unsigned int j=0; j<getInputDim(); j++)
      {
	for(unsigned int i=0; i<getNumData(); i++)
      	{
          g.setVal(gXorW.getVal(i, j), 0, counter);
          counter++;
        }
      }
    }
    else
    {
      for(unsigned int i=0; i<pbackModel->getOptNumParams(); i++)
      {
        g.setVal(gXorW.getVal(0, i), 0, counter);
      }
      throw ndlexceptions::NotImplementedError("Back constraints not yet implemented.");
    }
  }
  if(isSparseApproximation()) 
  {
    if(!isInducingFixed()) 
    {
      for(unsigned int j=0; j<getInputDim(); j++) 
      {
        for(unsigned int i=0; i<numActive; i++) 
        {
          g.setVal(gX_u.getVal(i, j), 0, counter);
          counter++;
        }
      }
    }
  }
  for(unsigned int i=0; i<pkern->getNumParams(); i++) 
  {
    g.setVal(g_param.getVal(0, i), 0, counter);
    counter++;
  }
  if(isOutputScaleLearnt()) 
  {
    for(unsigned int i=0; i<getOutputDim(); i++) 
    {
      g.setVal(g_scaleBias.getVal(0, i), 0, counter);
    }
    counter++;
  }
  if(isSparseApproximation()) 
  {
    double gBetaVal = gBeta.getVal(0);
    g.setVal(gBetaVal*betaTransform->gradfact(getBetaVal(0,0)), 0, counter);
    counter++;
  }
  //setADupToDate(false);
  return logLikelihood();

}
void CGp::updateG() const 
{
  if(!isMupToDate())
    throw ndlexceptions::Error("updateG() called when M is not updated.");
  unsigned int numKernParams = pkern->getNumParams();
  unsigned int numParams = numKernParams;
  
  updateK();
  updateAD();
  CMatrix tempG(1, numKernParams);
  CMatrix tempG2(1, numKernParams);
  CMatrix tmpV(getOutputDim(), 1);
  g_param.zeros();

  switch(getApproximationType()) 
  {
  case FTC:
    if(isOptimiseX()) // test if GP-LVM is used.
    {
      gXorW.zeros();
      pkern->getDiagGradX(dgKX, *pX);
    }
    for(unsigned int j=0; j<getOutputDim(); j++) 
    {
      updateCovGradient(j, tmpV); //covGrad = -(invK Y(:,j) Y(:,j)^t invK - invK)/2
      if(j==0) 
      {
        pkern->getGradTransParams(tempG, *pX, covGrad, true);
      }
      else 
      {
        pkern->getGradTransParams(tempG, *pX, covGrad, false);
      }
      for(unsigned int i=0; i<numKernParams; i++) 
      {
        g_param.addVal(tempG.getVal(i), i);
      }

      if(isOptimiseX())
      {
	for(unsigned int i=0; i<getNumData(); i++)
	{
	  pkern->getGradX(gKX, *pX, i, *pX); // gX[1...ndata].val(1:ndata,1:latentdim)
	  gKX.scale(2.0); // accounts for symmetric covariance
	  for(unsigned int l=0; l<getInputDim(); l++)
	  {
	    // deal with diagonal 
	    gKX.setVal(dgKX.getVal(i, l), i, l);
	  }
	  if(!isBackConstrained())
	  {
	    for(unsigned int k=0; k<getInputDim(); k++)
	    {
	      int ind = i + getNumData()*k;
	      gXorW.addVal(gKX.dotColCol(k, covGrad, i), ind);
	    }
	  }
	  else
	  {
	    throw ndlexceptions::NotImplementedError("Back constraints not yet implemented.");
	  }
	}
      }
    }
    break;
    
  case DTC:
  case DTCVAR:
  case FITC:
  case PITC:
    gpCovGrads();
    pkern->getGradTransParams(tempG, X_u, gK_uu, true);
    // should we regularize here? I.e. should the true at the end be false?
    pkern->getGradTransParams(tempG2, X_u, *pX, gK_uf, true);
    for(unsigned int i=0; i<numKernParams; i++) 
    {
      g_param.addVal(tempG.getVal(i), i);
      g_param.addVal(tempG2.getVal(i), i);
    }
    if(!isInducingFixed()) 
    {
      // Compute gradients with respect to X_u
      pkern->getDiagGradX(dgKX, X_u);
      // compute portion associated with gK_uf
      for(unsigned int i=0; i<numActive; i++) 
      {
        pkern->getGradX(gKX_uf, X_u, i, *pX);
        pkern->getGradX(gKX, X_u, i, X_u);
        gKX.scale(2.0); // accounts for fact that covGrad is symmetric.
        gKX.copyRowRow(i, dgKX, i);
        gX_u.resize(numActive, getInputDim());
        for(unsigned int j=0; j<getInputDim(); j++) 
        {
          gX_u.setVal(gKX.dotColCol(j, gK_uu, i), i, j);
        }
        for(unsigned int j=0; j<getInputDim(); j++) 
        {
          gX_u.addVal(gKX_uf.dotColRow(j, gK_uf, i), i, j); 
        }
      }
      if(isOptimiseX())
      {
	// for FITC and PITC need extra terms (see matlab code). This
	// is sufficient for DTC and DTCVAR.
	for(unsigned int i=0; i<getNumData(); i++)
	{
	  pkern->getGradX(gKX_uf2, *pX, i, X_u);
	  for(unsigned int j=0; j<getInputDim(); j++)
	  {
	    gXorW.setVal(gKX_uf2.dotColCol(j, gK_uf, i), i, j);
	  }
	}
      }
    }
    break;
  }
  switch(getApproximationType()) 
  {
  case FTC:
    break;    
  case DTC:
    break;
  case FITC:
  case DTCVAR:
    // deal with diagonal term's affect on kernel parameters.
    // false here stops gradients of priors being added ...
    pkern->getDiagGradTransParams(tempG, *pX, gLambda, false);
    for(unsigned int i=0; i<numKernParams; i++) 
    {
      g_param.addVal(tempG.getVal(i), i);
    }	
    if(isOptimiseX())
    {
      // using gKX_uf here as it is correct storeage size, really computing gKXdiag though!!
      pkern->getDiagGradX(gKX_uf, *pX);
      for(unsigned int i=0; i<getNumData(); i++)
      {
	for(unsigned int j=0; j<getInputDim(); j++)
	{
	  gXorW.addVal(gLambda.getVal(i)*gKX_uf.getVal(i, j), i, j);
	}
      }
    }
    break;
  case PITC:
    if(isOptimiseX())
      throw ndlexceptions::NotImplementedError("Optimisation of X not yet implemnted for PTIC sparse approximations");

    break;
  }
  
  if(isOutputScaleLearnt()) 
  {
    CMatrix &invKm = tmpV;  tmpV.resize(m.getRows(), 1);
    //  Need to fix this for FITC DTC etc.
    for(unsigned int j=0; j<getOutputDim(); j++) 
    {
      // recomputing this again is inefficient (already done in the likelihood).
      invKm.symvColCol(0, invK, m, j, 1.0, 0.0, "u");      
      g_scaleBias.setVal(1/scale.getVal(j)*(invKm.dotColCol(0, m, j)-1), j);
    }
  }
}

void CGp::gpCovGrads() const
{
  switch(getApproximationType()) 
  {
  case FTC:
    throw ndlexceptions::Error("gpCovGrads should not be called for FTC models.");
    break;
    
  case DTC:
  case DTCVAR:
    if(isSpherical())
    {
      E.gemm(K_uf, m, 1.0, 0.0, "n", "n");
      
      EET.gemm(E, E, 1.0, 0.0, "n", "t");
      EET.setSymmetric(true);
      
      AinvEET.symm(Ainv, EET, 1.0, 0.0, "u", "l");
      
      AinvEETAinv.gemm(AinvEET, Ainv, 1.0, 0.0, "n", "n");
      AinvEETAinv.setSymmetric(true);
      
      gK_uu.deepCopy(invK_uu);
      gK_uu.axpy(Ainv,-1.0/getBetaVal());
      gK_uu.scale((double)getOutputDim());
      gK_uu.axpy(AinvEETAinv, -1.0);
      gK_uu.setSymmetric(true);
      if(getApproximationType()==DTCVAR)
      {
	invK_uuK_uf.gemm(invK_uu, K_uf, 1.0, 0.0, "n", "n");
	gK_uu.syrk(invK_uuK_uf, -getBetaVal()*(double)getOutputDim(), 1.0, "l", "n");
      }
      gK_uu.scale(0.5);
      AinvK_uf.symm(Ainv, K_uf, 1.0, 0.0, "l", "l");
      
      EMT.gemm(E, m, 1.0, 0.0, "n", "t");
      
      AinvEMT.symm(Ainv, EMT, 1.0, 0.0, "l", "l");
      
      gK_uf.gemm(AinvEET, AinvK_uf, 1.0, 0.0, "n", "n");
      gK_uf.axpy(AinvEMT, -1.0);
      gK_uf.scale(getBetaVal());
      gK_uf.negate();
      gK_uf.axpy(AinvK_uf, -(double)getOutputDim());
      if(getApproximationType()==DTCVAR)
      {
	gK_uf.axpy(invK_uuK_uf, getBetaVal()*(double)getOutputDim());
      }
      double gbetaVal = ((double)(getNumData()-numActive)/getBetaVal()); 
      double temp = 0.0;
      for(unsigned int i=0; i < Ainv.getRows(); i++)
	temp += Ainv.dotRowRow(i, K_uu, i);
      gbetaVal += temp/(getBetaVal()*getBetaVal());
      gbetaVal *= (double)getOutputDim();
      temp = 0.0;
      for(unsigned int i=0; i < AinvEETAinv.getRows(); i++)
	temp += AinvEETAinv.dotRowRow(i, K_uu, i);
      gbetaVal += temp/getBetaVal();
      //CMatrix StoreNK(getNumData(), getOutputDim());
      for(unsigned int j=0; j<getOutputDim(); j++)
	gbetaVal -= m.dotColCol(j, m, j);
      gbetaVal += AinvEET.trace();
      if(getApproximationType()==DTCVAR)
      {
	gbetaVal -= (double)getOutputDim()*diagD.sum()/getBetaVal();
      }
      gbetaVal *= 0.5;
      gBeta.setVal(gbetaVal, 0);
      if(getApproximationType()==DTCVAR)
      {
	gLambda.setVals(-0.5*(double)getOutputDim()*getBetaVal());
      }
    }
    else 
    {
      throw ndlexceptions::Error("Non-spherical implementations not yet in place for DTC and DTCVAR");
    }
    break;

  case FITC:
    if(isSpherical())
    {
      //double betaVal = beta.getVal(0, 0);
      // using V as temp storage here.
      V.deepCopy(K_uf);
      for(unsigned int i=0; i<getNumData(); i++)
      {
	V.scaleCol(i, 1/diagD.getVal(i));
      }
      E.gemm(V, m, 1.0, 0.0, "n", "n");
      
      EET.gemm(E, E, 1.0, 0.0, "n", "t");
      EET.setSymmetric(true);
      AinvE.symm(Ainv, E, 1.0, 0.0, "l", "l");
      EMT.gemm(E, m, 1.0, 0.0, "n", "t");
      AinvEMT.symm(Ainv, EMT, 1.0, 0.0, "l", "l");
      K_ufdotTimesAinvEMT.deepCopy(AinvEMT);
      K_ufdotTimesAinvEMT.dotMultiply(K_uf);
      diagK_ufAinvEMT.sumCol(K_ufdotTimesAinvEMT, 1.0, 0.0);
      
      AinvEETAinv.gemm(AinvE, AinvE, 1.0, 0.0, "n", "t");
      AinvEETAinv.setSymmetric(true);
      // using Am as temp storage here.
      Am.deepCopy(Ainv);
      Am.scale((double)getOutputDim());
      Am.axpy(AinvEETAinv, getBetaVal());
      Am.setSymmetric(true);
      // using V as temp storage here.
      V.symm(Am, K_uf, 1.0, 0.0, "l", "l");
      V.dotMultiply(K_uf);
      diagK_ufdAinvplusAinvEETAinvK_fu.sumCol(V, 1.0, 0.0);
      invK_uuK_uf.symm(invK_uu, K_uf, 1.0, 0.0, "l", "l");
      invK_uuK_ufDinv.deepCopy(invK_uuK_uf);
      for(unsigned int i=0; i<getNumData(); i++)
      {
	invK_uuK_ufDinv.scaleCol(i, 1/diagD.getVal(i));
      }
      for(unsigned int i=0; i<getNumData(); i++)
	diagMMT.setVal(m.dotRowRow(i, m, i), i); 
      diagQ.deepCopy(diagK_ufdAinvplusAinvEETAinvK_fu);
      diagQ.axpy(diagD, -(double)getOutputDim());
      diagQ.axpy(diagMMT, getBetaVal());
      diagQ.axpy(diagK_ufAinvEMT, -2.0*getBetaVal());
      invK_uuK_ufDinvQ.deepCopy(invK_uuK_ufDinv);
      for(unsigned int i=0; i<getNumData(); i++)
      {
	invK_uuK_ufDinvQ.scaleCol(i, diagQ.getVal(i));
      }
      gK_uu.deepCopy(invK_uu);
      
      gK_uu.axpy(Ainv, -1.0/getBetaVal());
      gK_uu.scale((double)getOutputDim());
      gK_uu.axpy(AinvEETAinv, -1.0);
      gK_uu.gemm(invK_uuK_ufDinvQ, invK_uuK_ufDinv, getBetaVal(), 1.0, "n", "t");
      gK_uu.scale(0.5);	 
      gK_uu.setSymmetric(true);
 
      gK_uf.deepCopy(invK_uuK_ufDinvQ);
      gK_uf.symm(Ainv, K_uf, -(double)getOutputDim(), -getBetaVal(), "l", "l");
      gK_uf.symm(AinvEETAinv, K_uf, -getBetaVal(), 1.0, "l", "l");
      gK_uf.axpy(AinvEMT, getBetaVal());
      for(unsigned int j=0; j<getNumData(); j++)
      {
	gK_uf.scaleCol(j, 1/diagD.getVal(j));
      }
      gLambda.deepCopy(diagQ);
      gLambda.dotDivide(diagD);
      gLambda.scale(0.5*getBetaVal());
      gLambda.dotDivide(diagD);
      gBeta.setVal(-gLambda.sum()/(getBetaVal()*getBetaVal()), 0, 0);

    }
    else 
    {
      throw ndlexceptions::Error("Non-spherical implementations not yet in place for FITC");
    }
    break;
    
  case PITC:
    if(isSpherical())
    {
    }
    else 
    {
      throw ndlexceptions::Error("Non-spherical implementations not yet in place for PITC");
    }
    break;
  }
  
}
void CGp::pointLogLikelihood(const CMatrix& y, const CMatrix& X) const
{
  
}
#ifdef _NDLMATLAB
CGp::CGp(CMatrix* pinData, 
	 CMatrix* ptargetData, 
	 CKern* pkernel, 
	 CNoise* pnoiseModel, 
	 const string gpInfoFile, 
	 const string gpInfoVariable, int verbos) : 
  CMapModel(pinData->getCols(), ptargetData->getCols(), ptargetData->getRows()), CProbabilisticOptimisable(), 
  pX(pinData), py(ptargetData), 
  pkern(pkernel), pnoise(pnoiseModel)
{
  _init();
  initKernStoreage();
  setVerbosity(verbos);
  readMatlabFile(gpInfoFile, gpInfoVariable);
  updateK();
  updateAD();
}
mxArray* CGp::toMxArray() const
{
  int dims[1];
  dims[0]=1;
  mxArray* matlabArray;
  if(isSparseApproximation()) 
  {
    if(isInducingFixed()) 
    {
      const char* fieldNames[]={"d", "q", "k", "N", "learnScales", "approx", "fixInducing", "inducingIndices", "scale", "bias", "beta", "betaTransform", "type"};
      matlabArray = mxCreateStructArray(1, dims, 13, fieldNames);
      mxSetField(matlabArray, 0, "learnScales", convertMxArray(isOutputScaleLearnt()));
      mxSetField(matlabArray, 0, "approx", convertMxArray(getApproximationStr()));
      mxSetField(matlabArray, 0, "fixInducing", convertMxArray(true));
      // The inducing indices.
      vector<int> inducingMatlab = inducingIndices;
      for(unsigned int i=0; i<inducingMatlab.size(); i++)
	inducingMatlab[i]++;
      mxSetField(matlabArray, 0, "inducingIndices", convertMxArray(inducingMatlab));
      mxSetField(matlabArray, 0, "scale", scale.toMxArray());
      mxSetField(matlabArray, 0, "bias", bias.toMxArray());
      mxSetField(matlabArray, 0, "beta", beta.toMxArray());
      mxSetField(matlabArray, 0, "betaTransform", convertMxArray(betaTransform->getType()));
    }
    else 
    {
      const char* fieldNames[]={"d", "q", "k", "N", "learnScales", "approx", "fixInducing", "X_u", "scale", "bias", "beta", "betaTransform", "type"};
      matlabArray = mxCreateStructArray(1, dims, 13, fieldNames);
      
      mxSetField(matlabArray, 0, "learnScales", convertMxArray(isOutputScaleLearnt()));
      mxSetField(matlabArray, 0, "approx", convertMxArray(getApproximationStr()));
      mxSetField(matlabArray, 0, "fixInducing", convertMxArray(false));
      mxSetField(matlabArray, 0, "X_u", X_u.toMxArray());
      mxSetField(matlabArray, 0, "scale", scale.toMxArray());
      mxSetField(matlabArray, 0, "bias", bias.toMxArray());
      mxSetField(matlabArray, 0, "beta", beta.toMxArray());
      mxSetField(matlabArray, 0, "betaTransform", convertMxArray(betaTransform->getType()));
    }  
  }
  else {
    const char* fieldNames[]={"d", "q", "k", "N", "learnScales", "approx", "scale", "bias", "type"};
    matlabArray = mxCreateStructArray(1, dims, 9, fieldNames);
    mxSetField(matlabArray, 0, "learnScales", convertMxArray(isOutputScaleLearnt()));
    mxSetField(matlabArray, 0, "approx", convertMxArray(getApproximationStr()));
    mxSetField(matlabArray, 0, "scale", scale.toMxArray());
    mxSetField(matlabArray, 0, "bias", bias.toMxArray());
  }
  mxSetField(matlabArray, 0, "d", convertMxArray((double)getOutputDim()));
  mxSetField(matlabArray, 0, "q", convertMxArray((double)getInputDim()));
  mxSetField(matlabArray, 0, "k", convertMxArray((double)numActive));
  mxSetField(matlabArray, 0, "N", convertMxArray((double)getNumData()));
  mxSetField(matlabArray, 0, "type", convertMxArray(getType()));
  return matlabArray;
}
void CGp::fromMxArray(const mxArray* matlabArray)
{
  string mxType = mxArrayExtractStringField(matlabArray, "type");
  if(mxType!=getType())
  {
    throw ndlexceptions::FileReadError("Error mismatch between saved type, " + mxType + ", and Class type, " + getType() + ".");
  }
  setOutputScaleLearnt(mxArrayExtractBoolField(matlabArray, "learnScales"));
  setApproximationStr(mxArrayExtractStringField(matlabArray, "approx"));
  setOutputDim(mxArrayExtractIntField(matlabArray, "d"));
  setInputDim(mxArrayExtractIntField(matlabArray, "q"));
  setNumData(mxArrayExtractIntField(matlabArray, "N"));
  numActive = mxArrayExtractIntField(matlabArray, "k");
  initStoreage();
  if(isSparseApproximation()) 
  {
    setInducingFixed(mxArrayExtractBoolField(matlabArray, "fixInducing"));
    if(isInducingFixed()) 
    {
      inducingIndices = mxArrayExtractVectorIntField(matlabArray, "inducingIndices");
      numActive = inducingIndices.size();
      for(unsigned int i=0; i<inducingIndices.size(); i++)
	inducingIndices[i]--;
    }
    else 
    {
      X_u.fromMxArray(mxArrayExtractMxArrayField(matlabArray, "X_u"));
      DIMENSIONMATCH(numActive==X_u.getRows());
      DIMENSIONMATCH(getInputDim()==X_u.getCols());
    }
    beta.fromMxArray(mxArrayExtractMxArrayField(matlabArray, "beta"));
    string betaTransformType = mxArrayExtractStringField(matlabArray, "betaTransform");
    if(betaTransformType == "exp")
      betaTransform = new CExpTransform();
    else if(betaTransformType == "negLogLogit")
      betaTransform = new CNegLogLogitTransform();
    else
      ndlexceptions::NotImplementedError("fromMxArray: Not yet implemented transform types apart from exp and negLogLogit."); 
  }
  
  scale.fromMxArray(mxArrayExtractMxArrayField(matlabArray, "scale"));
  bias.fromMxArray(mxArrayExtractMxArrayField(matlabArray, "bias"));
  updateM();
}
#else /* not _NDLMATLAB */
#endif
// Optimise the GP with respect to latent positions and kernel parameters.
void CGp::optimise(unsigned int iters)
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

bool CGp::equals(const CGp& model, const double tol) const 
{
  if(!pnoise->equals(*model.pnoise, tol))
    return false;
  if(!pkern->equals(*model.pkern, tol))
    return false;
  if(!m.equals(model.m, tol))
    return false;
  if(isSparseApproximation())
    if(abs(getBetaVal()-model.getBetaVal())>tol)
      return false;
  if(!bias.equals(model.bias, tol))
    return false;
  if(!scale.equals(model.scale, tol))
    return false;
  if(!pX->equals(*model.pX, tol))
    return false;
  if(approximationType!=model.getApproximationType())
    return false;
  if(isSparseApproximation()) {
    if(numActive!=model.getNumActive())
      return false;
    if(!X_u.equals(model.X_u, tol))
      return false;
  }
  return true;
}

void CGp::display(ostream& os) const 
{
  if(isSparseApproximation())
    cout << "Sparse Approximation GP Model:" << endl << "Approx type: " << getApproximationStr() <<endl;
  else
    cout << "Standard GP Model: " << endl;
  cout << "Optimiser: " << getDefaultOptimiserStr() << endl;
  cout << "Data Set Size: " << getNumData() << endl;
  cout << "Kernel Type: " << endl;
  cout << "Scales learnt: " << isOutputScaleLearnt() << endl;
  cout << "X learnt: " << isOptimiseX() << endl;
  cout << "Bias: " << bias << endl;
  cout << "Scale: " << scale << endl;
  pnoise->display(os);
  pkern->display(os);
  if(isSparseApproximation()) {
    cout << "Inducing fixed: " << isInducingFixed() << endl;
    cout << "Beta Value: " << getBetaVal() << endl;
  }
  if(py && pX)
    cout << "Log likelihood: " << logLikelihood() << endl;
}

void CGp::readParamsFromStream(istream& in)
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
  setApproximationType(readIntFromStream(in, "sparseApproximation"));
  numActive = readIntFromStream(in, "numActive");
  // have all the info for resizing storage.
  initStoreage();

  if(isSparseApproximation()) 
  {
    betaTransform = CTransform::defaultPositive();
    beta.fromStream(in);
  }


  setOutputScaleLearnt(readBoolFromStream(in, "learnScale"));
  setOutputBiasLearnt(readBoolFromStream(in, "learnBias"));

  // next lines are scale and bias
  scale.fromStream(in);
  bias.fromStream(in);

  // load kernel
  pkern = readKernFromStream(in);
  initKernStoreage();
  // load noise
  pnoise = readNoiseFromStream(in);


  if(isSparseApproximation()) 
  {
    setInducingFixed(readBoolFromStream(in, "fixInducing"));
    X_u.fromStream(in);
    if(X_u.getCols()!=getInputDim())
      throw ndlexceptions::StreamFormatError("inputDim", "X_u columns doesn't match input dimension.");
  }
  setMupToDate(false);
  // Don't update M as you need py to be set for that.
  //updateM();
}


void CGp::writeParamsToStream(ostream& out) const 
{
  writeToStream(out, "baseType", getBaseType());
  writeToStream(out, "type", getType());
  writeToStream(out, "numData", getNumData());
  writeToStream(out, "outputDim", getOutputDim());
  writeToStream(out, "inputDim", getInputDim());

  writeToStream(out, "sparseApproximation", getApproximationType());
  writeToStream(out, "numActive", getNumActive());
  if(isSparseApproximation())
    beta.toStream(out);
  writeToStream(out, "learnScale", isOutputScaleLearnt());
  writeToStream(out, "learnBias", isOutputBiasLearnt());
  scale.toStream(out);
  bias.toStream(out);
  pkern->toStream(out);
  pnoise->toStream(out);
  if(isSparseApproximation()) 
  {
    writeToStream(out, "fixInducing", isInducingFixed());
    X_u.toStream(out);
  }
  //py->toStream(out);
  //pX->toStream(out);

}
// Functions which operate on the object
void writeGpToStream(const CGp& model, ostream& out) 
{
  model.toStream(out);
}

void writeGpToFile(const CGp& model, const string modelFileName, const string comment)
{
  model.toFile(modelFileName, comment);
}

CGp* readGpFromStream(istream& in)
{
  CGp* pmodel = new CGp();
  pmodel->fromStream(in);
  return pmodel;
}
    
CGp* readGpFromFile(const string modelFileName, int verbosity)
{
  // File is m, beta, X
  if(verbosity>0)
    cout << "Loading model file." << endl;
  ifstream in(modelFileName.c_str());
  if(!in.is_open()) throw ndlexceptions::FileReadError(modelFileName);
  CGp* pmodel;
  try 
  {
    pmodel = readGpFromStream(in);
  }
  catch(ndlexceptions::StreamFormatError err) 
  {
    throw ndlexceptions::FileFormatError(modelFileName, err);
  }
  if(verbosity>0)
    cout << "... done." << endl;
  in.close();
  pmodel->setVerbosity(verbosity);
  return pmodel; 
}
