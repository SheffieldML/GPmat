#ifndef CNOISE_H
#define CNOISE_H
#include <iostream>
#include <cmath>
#include <climits>
#include <vector>
#include "ndlexceptions.h"
#include "ndlutil.h"
#include "ndlstrutil.h"
#include "CMatrix.h"
#include "CTransform.h"
#include "COptimisable.h"
#include "CDist.h"
#include "CKern.h"
using namespace std;

const string NOISEVERSION="0.1";
// The basic noise class. The class allows its parameters to be transformed or optimised. Also noise models can be loaded from matlab files.
class CNoise : public CTransformable, public COptimisable, public CStreamInterface, public CMatInterface {
 public:
  // constructors
  CNoise() {  }
  CNoise(CMatrix* pyin) : py(pyin)
  {
    setNumData(py->getRows());
    setOutputDim(py->getCols());
  }
  virtual ~CNoise(){}
  
  // pure virtual functions
  virtual void initStoreage()=0;
  virtual void initNames()=0;
  virtual void initVals()=0;
  virtual void initParams()=0;
  virtual ostream& display(ostream& os)=0;
  virtual void setParams(const CMatrix& X)=0;
  virtual void getParams(CMatrix& X) const=0;
  virtual void getGradParams(CMatrix& g) const=0;
  virtual void getGradInputs(double& dlnZ_dmu, double& dlnZ_dvs, unsigned int i, unsigned int j) const=0;
  virtual void getGradInputs(CMatrix& dlnZ_dmu, CMatrix& dlnZ_dvs) const
  {
    // gradient with respect to mu and varSigma of log likelihood.  
    DIMENSIONMATCH(dlnZ_dmu.dimensionsMatch(dlnZ_dvs));
    DIMENSIONMATCH(dlnZ_dmu.getRows()==nData);
    DIMENSIONMATCH(dlnZ_dmu.getCols()==nProcesses);
    
    double gmu;
    double gvs;
    for(unsigned int i=0; i<dlnZ_dmu.getRows(); i++)
    {
      for(unsigned int j=0; j<dlnZ_dmu.getCols(); j++)
      {
	getGradInputs(gmu, gvs, i, j);
	dlnZ_dmu.setVal(gmu, i, j);
	dlnZ_dvs.setVal(gvs, i, j);
      }
    }
  }
  // Nu and G are combinations of gradients with respect to the mean and variance of the noise model.
  virtual void getNuG(CMatrix& g, CMatrix& nu, unsigned int index) const;
  virtual void updateSites(CMatrix& m, CMatrix& beta, unsigned int actIndex, const CMatrix& g, const CMatrix& nu, unsigned int index) const;
  virtual void test(const CMatrix& muout, const CMatrix& varSigmaOut, const CMatrix& yTest) const=0;
  virtual void out(CMatrix& yPred, const CMatrix& muTest, const CMatrix& varSigmaTest) const=0;
  virtual void out(CMatrix& yPred, CMatrix& probOut, const CMatrix& muTest, const CMatrix& varSigmaTest) const=0;
  virtual void likelihoods(CMatrix& L, const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const=0;
  virtual double logLikelihood(const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const=0;
  virtual double logLikelihood() const=0;
  virtual void writeParamsToStream(ostream& out) const;
  virtual void readParamsFromStream(istream& in);
#ifdef _NDLMATLAB
  virtual mxArray* toMxArray() const;
  virtual void fromMxArray(const mxArray* matlabArray);
  // Adds parameters to the mxArray.
  virtual void addParamToMxArray(mxArray* matlabArray) const;
  // Gets the parameters from the mxArray.
  virtual void extractParamFromMxArray(const mxArray* matlabArray);
#endif
  // non virtual functions
  inline void setType(const string name)
  {
    type = name;
  }
  inline string getType() const
  {
    return type;
  }
  string getBaseType() const
  {
    return "noise";
  }
  inline string getName() const
  {
    return noiseName;
  }
  inline void setName(const string name)
  {
    noiseName = name;
  }
  
  //void getGradX(CMatrix& gX, const CMatrix& dmu, const CMatrix& cvs);
  
  // Functions arising from optimisable.
  unsigned int getOptNumParams() const
  {
    return getNumParams();
  }
  void getOptParams(CMatrix& params) const
  {
    getTransParams(params);
  }
  void setOptParams(const CMatrix& params)
  {
    setTransParams(params);
  }
  double computeObjectiveGradParams(CMatrix& g) const
  {
    getGradTransParams(g);
    g.negate();
    return -logLikelihood();
  }
  double computeObjectiveVal() const
  {
    return -logLikelihood();
  }
  
  
  inline bool isLogConcave() const
  {
    return logConcave;
  }
  inline bool isSpherical() const
  {
    return spherical;
  }
  inline bool isMissing() const
  {
    return missing;
  }
  inline unsigned int getNumParams() const
  {
    return nParams;
  }
  inline unsigned int getOutputDim() const
  {
    return nProcesses;
  }
  inline unsigned int getNumData() const
  {
    return nData;
  }
  virtual string getParamName(unsigned int index) const
  {
    BOUNDCHECK(index>=0);
    BOUNDCHECK(index<paramNames.size());
    return paramNames[index];
  }
  void setParamName(const string paramName, unsigned int index)
  {
    BOUNDCHECK(index>=0);
    BOUNDCHECK(index<getNumParams());
    if(paramNames.size() == index)
      paramNames.push_back(paramName);
    else 
    {
      if(paramNames.size()<index)
	paramNames.resize(index+1, "no name");
      paramNames[index] = paramName;
    }
  }
  
  virtual void setMu(double val, unsigned int i, unsigned int j)=0;
  virtual double getMu(unsigned int i, unsigned int j) const=0;
  virtual void setVarSigma(double val, unsigned int i, unsigned int j)=0;
  virtual double getVarSigma(unsigned int i, unsigned int j) const=0;
  virtual double getTarget(unsigned int i, unsigned int j) const=0;
  virtual void setTarget(CMatrix* vals)
  {
    py=vals;
    setNumData(vals->getRows());
    setOutputDim(vals->getCols());
    initStoreage();
    initNames();      
  }
  
  virtual void setMus(double val)
  {
    for(unsigned int i=0; i<getNumData(); i++)
      for(unsigned int j=0; j<getOutputDim(); j++)
	setMu(val, i, j);
  }
  virtual void setVarSigmas(double val)
  {
    for(unsigned int i=0; i<getNumData(); i++)
      for(unsigned int j=0; j<getOutputDim(); j++)
	setVarSigma(val, i, j);
  }
  
  virtual void setMus(const CMatrix& vals)
  {
    DIMENSIONMATCH(vals.getRows()==getNumData());
    DIMENSIONMATCH(vals.getCols()==getOutputDim());
    for(unsigned int i=0; i<getNumData(); i++)
      for(unsigned int j=0; j<getOutputDim(); j++)
	setMu(vals.getVal(i, j), i, j);
  }
  virtual void setVarSigmas(const CMatrix& vals)
  {
    DIMENSIONMATCH(vals.getRows()==getNumData());
    DIMENSIONMATCH(vals.getCols()==getOutputDim());
    for(unsigned int i=0; i<getNumData(); i++)
      for(unsigned int j=0; j<getOutputDim(); j++)
	setVarSigma(vals.getVal(i, j), i, j);
  }
  bool equals(const CNoise& noise, double tol=ndlutil::MATCHTOL) const;
#ifdef _NDLMATLAB
  virtual void setVarSigmas(const mxArray* matlabArray)
  {
    CMatrix varsig;
    varsig.fromMxArray(matlabArray);
    setVarSigmas(varsig);
  }
  virtual void setMus(const mxArray* matlabArray)
  {
    CMatrix mus;
    mus.fromMxArray(matlabArray);
    setMus(mus);
  }
  virtual mxArray* varSigmaToMxArray() const
  {
    return varSigma.toMxArray();
  }
  virtual mxArray* muToMxArray() const
  {
    return mu.toMxArray();
  }
  
  virtual mxArray* targetToMxArray() const
  {
    return py->toMxArray();
  }
#endif

  CMatrix* py;
  
 protected:
  CMatrix mu;
  CMatrix varSigma;
  
  
  inline void setLogConcave(const bool val)
  {
    logConcave = val;
  }
  inline void setMissing(const bool val) 
  {
    missing = val;
  }
  inline void setSpherical(const bool val)
  {
    spherical = val;
  }
  inline void setNumParams(unsigned int num) 
  {
    nParams = num;
  }
  inline void setOutputDim(unsigned int num) 
  {
    nProcesses = num;
  }
  
 protected:
  inline void setNumData(unsigned int num) 
  {
    nData = num;
  }

 private:
   void _init();
  bool spherical;
  bool logConcave;
  bool missing;
  vector<string> paramNames;
  unsigned int nParams;
  unsigned int nProcesses;
  unsigned int nData;
  string type;
  string noiseName;
};
// The Gaussian noise model as commonly used in regression.
class CGaussianNoise : public CNoise {
 public:  
  // constructors
  CGaussianNoise()
  {
    _init();
  }
  CGaussianNoise(CMatrix* pyin) 
  {
    _init();
    setTarget(pyin);
    initNames();
    initVals();
    initParams();
  }
  ~CGaussianNoise();
  
  void initStoreage();
  void initNames();
  void initVals();
  void initParams();
  
  ostream& display(ostream& os);
  void setParams(const CMatrix& params);
  void setParam(double val, unsigned int index);
  void getParams(CMatrix& params) const;
  double getParam(unsigned int index) const;
  void getGradParams(CMatrix& g) const;
  void getGradInputs(double& gmu, double& gvs, unsigned int i, unsigned int j) const;
  void getNuG(CMatrix& g, CMatrix& nu, unsigned int index) const;
  void updateSites(CMatrix& m, CMatrix& beta, unsigned int actIndex, const CMatrix& g, const CMatrix& nu, unsigned int index) const;
  void test(const CMatrix& muout, const CMatrix& varSigmaOut, const CMatrix& yTest) const;
  void out(CMatrix& yPred, const CMatrix& muTest, const CMatrix& varSigmaTest) const;
  void out(CMatrix& yPred, CMatrix& errorBarOut, const CMatrix& muTest, const CMatrix& varSigmaTest) const;
  void likelihoods(CMatrix& L, const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  double logLikelihood(const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  double logLikelihood() const 
  {
    return logLikelihood(mu, varSigma, *py);
  }
  //void  getGradX(CMatrix& gX, const CMatrix& dmu, const CMatrix& cvs);
  
  inline double getMu(unsigned int i, unsigned int j) const 
  {
    return mu.getVal(i, j);
  }
  inline void setMu(double val, unsigned int i, unsigned int j) 
  {
    mu.setVal(val, i, j);
  }
  inline double getVarSigma(unsigned int i, unsigned int j) const 
  {
    return varSigma.getVal(i, j);
  }
  inline void setVarSigma(double val, unsigned int i, unsigned int j) 
  {
    SANITYCHECK(!isnan(val));
    SANITYCHECK(val>=0);
    varSigma.setVal(val, i, j);
  }
  inline double getTarget(unsigned int i, unsigned int j) const 
  {
    return py->getVal(i, j);
  }
  double getBiasVal(unsigned int index) const 
  {
    BOUNDCHECK(index<getOutputDim());
    BOUNDCHECK(index>=0);
    return bias.getVal(0, index);
  }
  void setBiasVal(double val, unsigned int index) 
  {
    BOUNDCHECK(index<getOutputDim());
    BOUNDCHECK(index>=0);
    bias.setVal(val, 0, index);
  }
  void setBias(const CMatrix& bia) 
  {
    DIMENSIONMATCH(bia.getRows()==1);
    DIMENSIONMATCH(bia.getCols()==getOutputDim());
    bias.deepCopy(bia);
  }
  
#ifdef _NDLMATLAB
  // Adds parameters to the mxArray.
  void addParamToMxArray(mxArray* matlabArray) const;
  // Gets the parameters from the mxArray.
  void extractParamFromMxArray(const mxArray* matlabArray);
#endif
 private:

  void _init();
  
  double sigma2;
  CMatrix bias;
};

// A scaled Gaussian noise model.
class CScaleNoise : public CNoise {
 public:  
  // constructors
  CScaleNoise()
  {
    _init();
  }
  CScaleNoise(CMatrix* pyin)
  {
    _init();
    setTarget(pyin);
    initVals();
    initParams(); 
  }
  ~CScaleNoise();
  
  double getScale(unsigned int index) const
  {
    BOUNDCHECK(index<getOutputDim()&&index>=0);
    return scale.getVal(index);
  }
  void setScale(double val, unsigned int index)
  {
    BOUNDCHECK(index<getOutputDim()&&index>=0);
    scale.setVal(val, index);
  }
  double getBias(unsigned int index) const
  {
    BOUNDCHECK(index<getOutputDim()&&index>=0);
    return bias.getVal(index);
  }
  void setBias(double val, unsigned int index)
  {
    BOUNDCHECK(index<getOutputDim()&&index>=0);
    bias.setVal(val, index);
  }
  void initStoreage();
  void initNames();
  void initVals();
  void initParams();
  
  ostream& display(ostream& os);
  void setParams(const CMatrix& params);
  void setParam(double val, unsigned int index);
  void getParams(CMatrix& params) const;
  double getParam(unsigned int index) const;
  void getGradParams(CMatrix& g) const;
  void getGradInputs(double& gmu, double& gvs, unsigned int i, unsigned int j) const;
  void getNuG(CMatrix& g, CMatrix& nu, unsigned int index) const;
  void updateSites(CMatrix& m, CMatrix& beta, unsigned int actIndex, const CMatrix& g, const CMatrix& nu, unsigned int index) const;
  void test(const CMatrix& muout, const CMatrix& varSigmaOut, const CMatrix& yTest) const;
  void out(CMatrix& yPred, const CMatrix& muTest, const CMatrix& varSigmaTest) const;
  void out(CMatrix& yPred, CMatrix& errorBarOut, const CMatrix& muTest, const CMatrix& varSigmaTest) const;
  void likelihoods(CMatrix& L, const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  double logLikelihood(const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  double logLikelihood() const
  {
    return logLikelihood(mu, varSigma, *py);
  }
  //void  getGradX(CMatrix& gX, const CMatrix& dmu, const CMatrix& cvs);
  
  inline double getMu(unsigned int i, unsigned int j) const
  {
    return mu.getVal(i, j);
  }
  inline void setMu(double val, unsigned int i, unsigned int j) 
  {
    mu.setVal(val, i, j);
  }
  inline double getVarSigma(unsigned int i, unsigned int j) const
  {
    return varSigma.getVal(i, j);
  }
  inline void setVarSigma(double val, unsigned int i, unsigned int j) 
  {
    SANITYCHECK(!isnan(val));
    SANITYCHECK(val>=0);
    varSigma.setVal(val, i, j);
  }
  inline double getTarget(unsigned int i, unsigned int j) const
  {
    return py->getVal(i, j);
  }
 
#ifdef _NDLMATLAB
  // Adds parameters to the mxArray.
  void addParamToMxArray(mxArray* matlabArray) const;
  // Gets the parameters from the mxArray.
  void extractParamFromMxArray(const mxArray* matlabArray);
#endif
 private:
  void _init();
  double sigma2;
  CMatrix bias;
  CMatrix scale;
};

// The probit noise model often used for classification.
class CProbitNoise : public CNoise 
{
 public:  
  // constructors
  CProbitNoise()
  {
    _init();
  }
  CProbitNoise(CMatrix* pyin) 
  {
    _init();
    setTarget(pyin);
    initVals();
    initParams(); 
  }
  ~CProbitNoise();
  
  void initStoreage();
  void initNames();
  void initVals();
  void initParams();
  
  ostream& display(ostream& os);
  void setParams(const CMatrix& params);
  void setParam(double val, unsigned int index);
  void getParams(CMatrix& params) const;
  double getParam(unsigned int index) const;
  void getGradParams(CMatrix& g) const;
  void getGradInputs(double& gmu, double& gvs, unsigned int i, unsigned int j) const;
  void test(const CMatrix& muout, const CMatrix& varSigmaOut, const CMatrix& yTest) const;
  void out(CMatrix& yPred, const CMatrix& muTest, const CMatrix& varSigmaTest) const;
  void out(CMatrix& yPred, CMatrix& probOut, const CMatrix& muTest, const CMatrix& varSigmaTest) const;
  void likelihoods(CMatrix& L, const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  double logLikelihood(const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  double logLikelihood() const
  {
    return logLikelihood(mu, varSigma, *py);
  }
  //void  getGradX(CMatrix& gX, const CMatrix& dmu, const CMatrix& cvs);
  
  inline double getMu(unsigned int i, unsigned int j) const
  {
    return mu.getVal(i, j);
  }
  inline void setMu(double val, unsigned int i, unsigned int j) 
  {
    mu.setVal(val, i, j);
  }
  inline double getVarSigma(unsigned int i, unsigned int j) const
  {
    return varSigma.getVal(i, j);
  }
  inline void setVarSigma(double val, unsigned int i, unsigned int j) 
  {
    SANITYCHECK(!isnan(val));
    SANITYCHECK(val>=0);
    varSigma.setVal(val, i, j);
  }
  inline double getTarget(unsigned int i, unsigned int j) const
  {
    return py->getVal(i, j);
  }
  
#ifdef _NDLMATLAB
  // Adds parameters to the mxArray.
  void addParamToMxArray(mxArray* matlabArray) const;
  // Gets the parameters from the mxArray.
  void extractParamFromMxArray(const mxArray* matlabArray);
#endif
  
 private:
  void _init();
  double sigma2;
  CMatrix bias;
};

// The null category noise model for semi-supervised learning (see Lawrence and Jordan in NIPS 2004).
class CNcnmNoise : public CNoise {
 public:  
  // constructors
  CNcnmNoise()
  {
    _init();
  }
  CNcnmNoise(CMatrix* pyin) 
  {
    _init();
    setTarget(pyin);
    initVals();
    initParams(); 
  }
  ~CNcnmNoise();
  
  void initStoreage();
  void initNames();
  void initVals();
  void initParams();

  ostream& display(ostream& os);
  void setParams(const CMatrix& params);
  void setParam(double val, unsigned int index);
  void getParams(CMatrix& params) const;
  double getParam(unsigned int index) const;
  void getGradParams(CMatrix& g) const;
  void getGradInputs(double& gmu, double& gvs, unsigned int i, unsigned int j) const;
  void test(const CMatrix& muout, const CMatrix& varSigmaOut, const CMatrix& yTest) const;
  void out(CMatrix& yPred, const CMatrix& muTest, const CMatrix& varSigmaTest) const;
  void out(CMatrix& yPred, CMatrix& probOut, const CMatrix& muTest, const CMatrix& varSigmaTest) const;
  void likelihoods(CMatrix& L, const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  double logLikelihood(const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  void writeParamsToStream(ostream& out) const;
  void readParamsFromStream(istream& in);
  double logLikelihood() const
  {
    return logLikelihood(mu, varSigma, *py);
  }
  //void  getGradX(CMatrix& gX, const CMatrix& dmu, const CMatrix& cvs);
  
  inline double getMu(unsigned int i, unsigned int j) const
  {
    return mu.getVal(i, j);
  }
  inline void setMu(double val, unsigned int i, unsigned int j) 
  {
    mu.setVal(val, i, j);
  }
  inline double getVarSigma(unsigned int i, unsigned int j) const
  {
    return varSigma.getVal(i, j);
  }
  inline void setVarSigma(double val, unsigned int i, unsigned int j) 
  {
    SANITYCHECK(!isnan(val));
    SANITYCHECK(val>=0);
    varSigma.setVal(val, i, j);
  }
  inline void setSplitGamma(const bool val)
  {
    splitGamma=val;
  }
  bool isSplitGamma() const
  {
    return splitGamma;
  }
  inline double getTarget(unsigned int i, unsigned int j) const
  {
    return py->getVal(i, j);
  }
  void setTarget(CMatrix* vals)
  {
    setSplitGamma(false);
    CNoise::setTarget(vals);
  }
#ifdef _NDLMATLAB
  // Adds parameters to the mxArray.
  void addParamToMxArray(mxArray* matlabArray) const;
  // Gets the parameters from the mxArray.
  void extractParamFromMxArray(const mxArray* matlabArray);
#endif
  
 private:
  void _init();
  double gamman;
  double gammap;
  bool splitGamma;
  double width;
  double sigma2;
  CMatrix bias;
};
class COrderedNoise : public CNoise {
 public:  
  // constructors
  COrderedNoise()
  {
    _init();
  }
  COrderedNoise(CMatrix* pyin) : numCats(1)
  {
    _init();
    setTarget(pyin);
    initVals();
    initParams(); 
  }
  COrderedNoise(CMatrix* pyin, unsigned int numCts) : numCats(numCts)
  {
    _init();
    setTarget(pyin);
    initVals();
    initParams(); 
  }
  ~COrderedNoise();
  
  unsigned int getNumCategories() const
  {
    return numCats;
  }
  void setNumCategories(unsigned int val)
  {
    SANITYCHECK(val>1);
    numCats=val;
    initStoreage();
  }
  
  void initStoreage();
  void initNames();
  void initVals();
  void initParams();

  ostream& display(ostream& os);
  void setParams(const CMatrix& params);
  void setParam(double val, unsigned int index);
  void getParams(CMatrix& params) const;
  double getParam(unsigned int index) const;
  void getGradParams(CMatrix& g) const;
  void getGradInputs(double& gmu, double& gvs, unsigned int i, unsigned int j) const;
  void test(const CMatrix& muout, const CMatrix& varSigmaOut, const CMatrix& yTest) const;
  void out(CMatrix& yPred, const CMatrix& muTest, const CMatrix& varSigmaTest) const;
  void out(CMatrix& yPred, CMatrix& probOut, const CMatrix& muTest, const CMatrix& varSigmaTest) const;
  void likelihoods(CMatrix& L, const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  double logLikelihood(const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  void writeParamsToStream(ostream& out) const;
  void readParamsFromStream(istream& in);
  double logLikelihood() const
  {
    return logLikelihood(mu, varSigma, *py);
  }
  //void  getGradX(CMatrix& gX, const CMatrix& dmu, const CMatrix& cvs);
  
  inline double getMu(unsigned int i, unsigned int j) const
  {
    return mu.getVal(i, j);
  }
  inline void setMu(double val, unsigned int i, unsigned int j) 
  {
    mu.setVal(val, i, j);
  }
  inline double getVarSigma(unsigned int i, unsigned int j) const
  {
    return varSigma.getVal(i, j);
  }
  inline void setVarSigma(double val, unsigned int i, unsigned int j) 
  {
    SANITYCHECK(!isnan(val));
    SANITYCHECK(val>=0);
    varSigma.setVal(val, i, j);
  }
  inline double getTarget(unsigned int i, unsigned int j) const
  {
    return py->getVal(i, j);
  }
#ifdef _NDLMATLAB
  // Adds parameters to the mxArray.
  void addParamToMxArray(mxArray* matlabArray) const;
  // Gets the parameters from the mxArray.
  void extractParamFromMxArray(const mxArray* matlabArray);
#endif
  
 private:
   void _init();
  unsigned int numCats;
  CMatrix widths;
  double sigma2;
  CMatrix bias;
  
  mutable CMatrix gwidth;
};


void writeNoiseToStream(const CNoise& noise, ostream& out);
CNoise* readNoiseFromStream(istream& in);


#endif
