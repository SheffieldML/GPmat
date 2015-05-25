#ifndef CGP_H
#define CGP_H
#include "CMltools.h"

using namespace std;

const string GPVERSION="0.1";

class CGp : public CMapModel, public CProbabilisticOptimisable, public CStreamInterface, public CMatInterface
{
public:
  enum 
  {
    FTC, 
    DTC, 
    FITC, 
    PITC,
    DTCVAR
  };

  CGp();

  // Constructor given a kernel
  //  CGp(CKern& kernel, CScaleNoise& nois, CMatrix& Xin, int verbos=2);
  // Constructor given a kernel and sparse approximation
  CGp(CKern* kernel, CNoise* nois, CMatrix* Xin, int approxType=FTC, unsigned int actSetSize=0, int verbos=2);

  CGp(unsigned int q, unsigned int d, CMatrix* Xin, CMatrix* yin, CKern* kernel, CNoise* nois, int approxType=FTC, unsigned int actSetSize=0, int verbos=2);

#ifdef _NDLMATLAB
  // Constructor using file containing gpInfo.
  CGp(CMatrix* inData, CMatrix* targetData, CKern* kernel, CNoise* noiseModel, const string gpInfoFile, const string gpInfoVariable, int verbos=2);
#endif

  // initialise storeage associated with optimizing inputs.
  void initOptimiseXStoreage();
  // initialise storeage associated with active set size.
  void initSparseStoreage();
  // Initialise storeage for anything dependent on number of kernel parameters.
  void initKernStoreage();
  // Initialise the storeage for the model.
  void initStoreage();
  // Set the initial values for the model.
  void initVals();
  
  // For MapModel interface.
  void out(CMatrix& yPred, const CMatrix& inData) const;
  void out(CMatrix& yPred, CMatrix& probPred, const CMatrix& inData) const;
  double outGradParams(CMatrix& g, const CMatrix& Xin, unsigned int pointNo, unsigned int outputNo) const;
  double outGradX(CMatrix& g, const CMatrix& Xin, unsigned int pointNo, unsigned int outputNo) const;

  // update alpha representation.
  void updateAlpha() const;
  void posteriorMeanVar(CMatrix& mu, CMatrix& varSigma, const CMatrix& X) const;
  void posteriorMean(CMatrix& mu, const CMatrix& X) const;
  // Gradient routines
  void updateCovGradient(unsigned int index, CMatrix &work_invK_Y) const;


  virtual void updateX();
  void updateM() const;
  // update K and dynK and all derived quantities if they are dirty.
  void updateK() const;
  // Update A and D representations
  void updateAD() const;
  // update the gradient matrices.
  void updateG() const;

  // compute the approximation to the log likelihood.
  virtual double logLikelihood() const;
  // compute the gradients of the approximation wrt parameters.
  virtual double logLikelihoodGradient(CMatrix& g) const;
  void gpCovGrads() const;
  virtual void pointLogLikelihood(const CMatrix& y, const CMatrix& X) const;
  void optimise(unsigned int iters=1000);
  bool equals(const CGp& model, double tol=ndlutil::MATCHTOL) const;
  void display(ostream& os) const;
  
  virtual unsigned int getOptNumParams() const;
  virtual void getOptParams(CMatrix& param) const;
  virtual void setOptParams(const CMatrix& param);

  virtual string getNoiseType() const;

#ifdef _NDLMATLAB
  mxArray* toMxArray() const;
  void fromMxArray(const mxArray* matlabArray);
#endif

  void readParamsFromStream(istream& in);
  void writeParamsToStream(ostream& out) const;

  inline unsigned int getNumActive() const 
  {
    return numActive;
  }
  void setNumActive(unsigned int val) 
  {
    numActive = val;
  }
  
  void setTargetVals(CMatrix& yvals) 
  {
    DIMENSIONMATCH(yvals.getCols()==getOutputDim());
    DIMENSIONMATCH(yvals.getRows()==getNumData());
    py=&yvals;
    pnoise->setTarget(py);
  }
  void setTargetVals(CMatrix* pyvals) 
  {
    DIMENSIONMATCH(pyvals->getCols()==getOutputDim());
    DIMENSIONMATCH(pyvals->getRows()==getNumData());
    py=pyvals;
    pnoise->setTarget(py);
  }
  void setInputVals(CMatrix& Xvals) 
  {
    DIMENSIONMATCH(Xvals.getCols()==getInputDim());
    DIMENSIONMATCH(Xvals.getRows()==getNumData());
    pX=&Xvals;
  }
  void setInputVals(CMatrix* pXvals) 
  {
    DIMENSIONMATCH(pXvals->getCols()==getInputDim());
    DIMENSIONMATCH(pXvals->getRows()==getNumData());
    pX=pXvals;
  }
  void setInducingVals(CMatrix& Xvals) 
  {
    SANITYCHECK(isSparseApproximation());
    DIMENSIONMATCH(Xvals.getCols()==getInputDim());
    DIMENSIONMATCH(Xvals.getRows()==numActive);
    X_u.deepCopy(Xvals);
  }
  // Flag which indicates whether scales are to be learnt.
  bool isOutputScaleLearnt() const 
  {
    return outputScaleLearnt;
  }
  void setOutputScaleLearnt(const bool val) {
    outputScaleLearnt=val;
  }
  bool isOutputBiasLearnt() const 
  {
    return outputBiasLearnt;
  }
  void setOutputBiasLearnt(const bool val) {
    outputBiasLearnt=val;
  }
  double getScaleVal(unsigned int index) const {
    BOUNDCHECK(index<getOutputDim());
    return scale.getVal(0, index);
  }
  void setScaleVal(double val, unsigned int index) {
    BOUNDCHECK(index<getOutputDim());
    scale.setVal(val, 0, index);
    setMupToDate(false);
  }
  void setScale(const CMatrix& scal) {
    DIMENSIONMATCH(scal.getRows()==1);
    DIMENSIONMATCH(scal.getCols()==getOutputDim());
    scale.deepCopy(scal);
    setMupToDate(false);
  } 
  void setActiveSet(const CMatrix& Xu) {
    DIMENSIONMATCH(Xu.getRows()==numActive);
    DIMENSIONMATCH(Xu.getCols()==getInputDim());
    X_u.deepCopy(Xu);
    setKupToDate(false);
  }
  double getBiasVal(unsigned int index) const {
    BOUNDCHECK(index<getOutputDim());
    return bias.getVal(0, index);
  }
  void setBiasVal(double val, unsigned int index) {
    BOUNDCHECK(index<getOutputDim());
    bias.setVal(val, 0, index);
    setMupToDate(false);
  }
  void setBias(const CMatrix& bia) {
    DIMENSIONMATCH(bia.getRows()==1);
    DIMENSIONMATCH(bia.getCols()==getOutputDim());
    bias.deepCopy(bia);
    setMupToDate(false);
  }
  double getBetaVal(unsigned int i=0, unsigned int j=0) const {
    BOUNDCHECK(i<beta.getRows());
    BOUNDCHECK(j<beta.getCols());
    return beta.getVal(i, j);
  }
  void setBetaVal(double val, unsigned int i=0, unsigned int j=0) {
    BOUNDCHECK(i<beta.getRows());
    BOUNDCHECK(j<beta.getCols());
    beta.setVal(val, i, j);
    setADupToDate(false);
  }
  void setBetaVals(double val) {
    beta.setVals(val);
    setADupToDate(false);
  }
  void setBeta(const CMatrix& bet) {
    DIMENSIONMATCH(bet.dimensionsMatch(beta));
    beta.deepCopy(bet);
    setADupToDate(false);
  }
  int getApproximationType() const {
    return approximationType;
  }
  string getApproximationStr() const {
    switch(approximationType) {
    case FTC:
      return "ftc";
    case DTC:
      return "dtc";
    case DTCVAR:
      return "dtcvar";
    case FITC:
      return "fitc";
    case PITC:
      return "pitc";
    default:
      throw ndlexceptions::Error("Unknown approximation type");
    }
  }
  void setApproximationStr(const string val) {
    if(val=="ftc") 
      setApproximationType(FTC);
    else if(val=="dtc")
      setApproximationType(DTC);
    else if(val=="dtcvar")
      setApproximationType(DTCVAR);
    else if(val=="fitc")
      setApproximationType(FITC);
    else if(val=="pitc")
      setApproximationType(PITC);
    else
      throw ndlexceptions::Error("Unknown approximation type");
  }

  void setApproximationType(unsigned int val) {
    approximationType=val;
    if(approximationType == FTC)
      setSparseApproximation(false);
    else if(approximationType == DTC)
      setSparseApproximation(true);
    else if(approximationType == DTCVAR)
      setSparseApproximation(true);
    else if(approximationType == FITC)
      setSparseApproximation(true);
    else if (approximationType == PITC)
      setSparseApproximation(true);
    else
      throw ndlexceptions::Error("Unknown approximation type");
  }
  // Flag which indicates if a sparse approximation is used.
  bool isSparseApproximation() const 
  {
    return sparseApproximation;
  }
  void setSparseApproximation(const bool val) 
  {
    sparseApproximation=val;
  }
  // Whether noise model on outputs is spherical.
  bool isSpherical() const 
  {
    return spherical;
  }
  void setSpherical(const bool val) {
    spherical = val;
  }
  // Whether inducing variables are fixed.
  bool isInducingFixed() const 
  {
    return inducingFixed;
  }
  void setInducingFixed(const bool val) 
  {
    inducingFixed = val;
  }
  // Flag which indicates if K/Kinv/DynK/DynKInv need recomputation.
  bool isKupToDate() const 
  {
    return KupToDate;
  }
  void setKupToDate(const bool val) const 
  {
    KupToDate = val;
    if(!KupToDate)
    {
      setAlphaUpToDate(false);
      setADupToDate(false);
    }
  }
  // Flag which indicates if A,D etc. need recomputation.
  bool isADupToDate() const 
  {
    return ADupToDate;
  }
  void setADupToDate(const bool val) const 
  {
    ADupToDate = val;
  }
  // Flag which indicates if alphas need recomputation.
  bool isAlphaUpToDate() const 
  {
    return AlphaUpToDate;
  }
  void setAlphaUpToDate(const bool val) const 
  {
    AlphaUpToDate = val;
  }
  // Flag which indicates if M needs recomputation.
  bool isMupToDate() const 
  {
    return MupToDate;
  }
  void setMupToDate(const bool val) const 
  {
    MupToDate = val;
    if(!MupToDate)
    {
      setAlphaUpToDate(false);
      setADupToDate(false);
    }
  }
  bool isOptimiseX() const 
  {
    return optimiseX;
  }
  void setOptimiseX(const bool val) 
  {
    bool change = false;
    if(val!=optimiseX)
      change = true;
    optimiseX = val;
    if(change)
      initOptimiseXStoreage();
  }
  bool isBackConstrained() const
  {
    return backConstrained;
  }
  void setBackConstrained(const bool val)
  {
    backConstrained = val;
  }
 
  CMapModel* backConstraintModel; // for mapping constraints on latent variables.
  CMatrix X_u; // for inducing variables if needed.

  CMatrix* pX;
  CMatrix* py;  // target data.

  // if debugging, make lots of these variables available for checking in python.
#ifndef DBG  
 private:
#endif
  bool optimiseX;
  bool backConstrained;

  mutable CMatrix m;  // scaled and biased Y
  mutable CMatrix Alpha; // SVM style 'alphas'.
  

  // Temporary variables for sparse approximations.
  mutable CMatrix Am;  
  mutable CMatrix Lm;  
  mutable CMatrix invLmV;  
  mutable CMatrix bet;  
  mutable CMatrix diagD;  
  mutable CMatrix sqrtDiagD;  
  mutable CMatrix scaledM; // numData*outputDim storage matrix.
  mutable CMatrix V;  // numActive*numData storage matrix
  mutable CMatrix E; // numActive*numData storage matrix
  mutable CMatrix EET; // numActive*numActive storage matrix.
  mutable CMatrix AinvE; // numActive*outputDim storage matrix.
  mutable CMatrix AinvEET; // numActive*numActive storage matrix.
  mutable CMatrix AinvEETAinv; // numActive*numActive storage matrix.
  mutable CMatrix AinvK_uf; // numActive*numData storage matrix.
  mutable CMatrix EMT; // numActive*numData storage matrix.
  mutable CMatrix AinvEMT; // numActive*numData storage matrix.
  mutable CMatrix invK_uuK_uf; // numActive*numData storage matrix.
  mutable CMatrix invK_uuK_ufDinv; // numActive*numData storage matrix. 
  mutable CMatrix invK_uuK_ufDinvQ; // numActive*numData storage matrix. 
  mutable CMatrix diagMMT; // numData*1 storage matrix.
  mutable CMatrix diagQ; // numData*1 storage matrix.
  mutable CMatrix diagK_ufdAinvplusAinvEETAinvK_fu; // numData*1 storage matrix.
  mutable CMatrix K_ufdotTimesAinvEMT; // numActive*numData storage matrix.
  mutable CMatrix diagK_ufAinvEMT; // numData*1 storage matrix.
  
  
  CMatrix beta;
  CMatrix nu;
  CMatrix scale;
  CMatrix bias;
  
  CMatrix g;
  CKern* pkern;
  CNoise* pnoise;
  
  unsigned int maxTries;

  mutable vector<int> inducingIndices;
  mutable CMatrix gDiagX;
  mutable CMatrix gK_uf;
  mutable CMatrix gK_uu;
  mutable CMatrix gK_star;
  mutable CMatrix gLambda;
  mutable CMatrix gBeta;

  mutable CMatrix K;
  mutable CMatrix invK;
  mutable CMatrix diagK;
  mutable CMatrix LcholK;
  mutable double logDetK;

  mutable CMatrix covGrad;
  mutable CMatrix tempgX;

  mutable CMatrix K_uu;
  mutable CMatrix invK_uu;
  mutable double logDetK_uu;

  mutable CMatrix StorekN;

  mutable CMatrix K_uf;

  mutable CMatrix A;
  mutable CMatrix Ainv;
  mutable CMatrix LcholA;
  mutable double logDetA;

  // Representations for intermediate gradients.
  mutable CMatrix dgKX;
  mutable CMatrix gKX;
  mutable CMatrix gKX_uf;
  mutable CMatrix gKX_uf2;

  // Matrices where gradients are temporarily stored.
  mutable CMatrix g_scaleBias;
  mutable CMatrix g_param; 
  mutable CMatrix gX_u;
  mutable CMatrix gXorW;
  CMapModel* pbackModel;
  // if debugging, still make the remainder available.
#ifdef DBG
 private:
#endif
  void _init();
  void _updateK() const; // update K with the inverse of the kernel plus beta terms computed from the active points.
  void _updateInvK(unsigned int dim=0) const;
  void _testComputeKx(CMatrix& kX, const CMatrix& Xin) const; /// compute the kernel for some test points.
  
  void _posteriorMean(CMatrix& mu, const CMatrix& kX) const; /// compute the posterior mean given the kernel evaluated at test points.
  void _posteriorVar(CMatrix& varSigma, CMatrix& kX, const CMatrix& Xin) const; /// compute the posterior variance given the kernel evaluated at test points.
  double jitter;

  CTransform* betaTransform;

  unsigned int numActive;
  int approximationType; /// FTC, DTC, FITC, PITC, DTCVAR
  int numCovStruct;

  bool outputScaleLearnt;
  bool outputBiasLearnt;
  bool sparseApproximation;
  bool terminate;
  bool epUpdate;
  bool loadedModel;
  bool spherical;
  bool inducingFixed;
  mutable bool KupToDate;
  mutable bool AlphaUpToDate;
  mutable bool ADupToDate;
  mutable bool MupToDate;
  
  string type;
};

// Functions which operate on the object
void writeGpToStream(const CGp& model, ostream& out);
void writeGpToFile(const CGp& model, const string modelFileName, const string comment="");
CGp* readGpFromStream(istream& in);
CGp* readGpFromFile(const string modelfileName, int verbosity=2);


#endif /* CGP_H */
