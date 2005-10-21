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
class CNoise : public CTransformable, public COptimisable, public CMatinterface {
 public:
  // constructors
  CNoise() {}
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
  virtual void getGradInputs(double& dlnZ_dmu, double& dlnZ_dvs, int i, int j) const=0;
  virtual void getGradInputs(CMatrix& dlnZ_dmu, CMatrix& dlnZ_dvs) const
    {
      // gradient with respect to mu and varSigma of log likelihood.  
      assert(dlnZ_dmu.dimensionsMatch(dlnZ_dvs));
      assert(dlnZ_dmu.getRows()==nData);
      assert(dlnZ_dmu.getCols()==nProcesses);
      
      double gmu;
      double gvs;
      for(int i=0; i<dlnZ_dmu.getRows(); i++)
	{
	  for(int j=0; j<dlnZ_dmu.getCols(); j++)
	    {
	      getGradInputs(gmu, gvs, i, j);
	      dlnZ_dmu.setVal(gmu, i, j);
	      dlnZ_dvs.setVal(gvs, i, j);
	    }
	}
    }
  // Nu and G are combinations of gradients with respect to the mean and variance of the noise model.
  virtual void getNuG(CMatrix& g, CMatrix& nu, int index) const;
  virtual void updateSites(CMatrix& m, CMatrix& beta, int actIndex, const CMatrix& g, const CMatrix& nu, int index) const;
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
  inline void setVerbosity(int val)
    {
      verbosity = val;
    }
  inline int getVerbosity() const
    {
      return verbosity;
    }
  inline void setType(const string name)
    {
      type = name;
    }
  inline string getType() const
    {
      return type;
    }
  inline string getNoiseName() const
    {
      return noiseName;
    }
  inline void setNoiseName(const string name)
    {
      noiseName = name;
    }

  void getGradX(CMatrix& gX, const CMatrix& dmu, const CMatrix& cvs);

  // Functions arising from optimisable.
  int getOptNumParams() const
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
  void computeObjectiveGradParams(CMatrix& g) const
    {
      getGradTransParams(g);
      g.negate();
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
  inline int getNumParams() const
    {
      return nParams;
    }
  inline int getNumProcesses() const
    {
      return nProcesses;
    }
  inline int getNumData() const
    {
      return nData;
    }
  virtual string getParamName(int index) const
    {
      assert(index>=0);
      assert(index<paramNames.size());
      return paramNames[index];
    }
  void setParamName(const string paramName, int index)
    {
      assert(index>=0);
      assert(index<getNumParams());
      if(paramNames.size() == index)
	paramNames.push_back(paramName);
      else 
	{
	  if(paramNames.size()<index)
	    paramNames.resize(index+1, "no name");
	  paramNames[index] = paramName;
	}
    }

  virtual void setMu(double val, int i, int j)=0;
  virtual double getMu(int i, int j) const=0;
  virtual void setVarSigma(double val, int i, int j)=0;
  virtual double getVarSigma(int i, int j) const=0;
  virtual double getTarget(int i, int j) const=0;
  virtual void setTarget(const CMatrix& vals)=0;

  virtual void setMus(double val)
    {
      for(int i=0; i<getNumData(); i++)
	for(int j=0; j<getNumProcesses(); j++)
	  setMu(val, i, j);
    }
  virtual void setVarSigmas(double val)
    {
      for(int i=0; i<getNumData(); i++)
	for(int j=0; j<getNumProcesses(); j++)
	  setVarSigma(val, i, j);
    }

  virtual void setMus(const CMatrix& vals)
    {
      assert(vals.getRows()==getNumData());
      assert(vals.getCols()==getNumProcesses());
      for(int i=0; i<getNumData(); i++)
	for(int j=0; j<getNumProcesses(); j++)
	  setMu(vals.getVal(i, j), i, j);
    }
  virtual void setVarSigmas(const CMatrix& vals)
    {
      assert(vals.getRows()==getNumData());
      assert(vals.getCols()==getNumProcesses());
      for(int i=0; i<getNumData(); i++)
	for(int j=0; j<getNumProcesses(); j++)
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
      return y.toMxArray();
    }
#endif
 
 protected:
  CMatrix mu;
  CMatrix varSigma;
  CMatrix y;


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
  inline void setNumParams(int num) 
    {
      nParams = num;
    }
  inline void setNumProcesses(int num) 
    {
      nProcesses = num;
      initStoreage();
      initNames();      
    }
  
 protected:
  inline void setNumData(int num) 
    {
      nData = num;
    }

 private:
  int verbosity;
  bool spherical;
  bool logConcave;
  bool missing;
  vector<string> paramNames;
  int nParams;
  int nProcesses;
  int nData;
  string type;
  string noiseName;
};
// The Gaussian noise model as commonly used in regression.
class CGaussianNoise : public CNoise {
 public:  
  // constructors
  CGaussianNoise(){}
  CGaussianNoise(const CMatrix& yin)
    {
      setTarget(yin);
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
  void setParam(double val, int index);
  void getParams(CMatrix& params) const;
  double getParam(int index) const;
  void getGradParams(CMatrix& g) const;
  void getGradInputs(double& gmu, double& gvs, int i, int j) const;
  void getNuG(CMatrix& g, CMatrix& nu, int index) const;
  void updateSites(CMatrix& m, CMatrix& beta, int actIndex, const CMatrix& g, const CMatrix& nu, int index) const;
  void test(const CMatrix& muout, const CMatrix& varSigmaOut, const CMatrix& yTest) const;
  void out(CMatrix& yPred, const CMatrix& muTest, const CMatrix& varSigmaTest) const;
  void out(CMatrix& yPred, CMatrix& errorBarOut, const CMatrix& muTest, const CMatrix& varSigmaTest) const;
  void likelihoods(CMatrix& L, const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  double logLikelihood(const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  double logLikelihood() const
    {
      return logLikelihood(mu, varSigma, y);
    }
  void  getGradX(CMatrix& gX, const CMatrix& dmu, const CMatrix& cvs);
  
  inline double getMu(int i, int j) const
    {
      return mu.getVal(i, j);
    }
  inline void setMu(double val, int i, int j) 
    {
      mu.setVal(val, i, j);
    }
  inline double getVarSigma(int i, int j) const
    {
      return varSigma.getVal(i, j);
    }
  inline void setVarSigma(double val, int i, int j) 
    {
      assert(!isnan(val));
      assert(val>=0);
      varSigma.setVal(val, i, j);
    }
  inline double getTarget(int i, int j) const
    {
      return y.getVal(i, j);
    }
  void setTarget(const CMatrix& vals)
    {
      y.deepCopy(vals);
      setNumData(vals.getRows());
      setNumProcesses(vals.getCols());
    }
#ifdef _NDLMATLAB
  // Adds parameters to the mxArray.
  void addParamToMxArray(mxArray* matlabArray) const;
  // Gets the parameters from the mxArray.
  void extractParamFromMxArray(const mxArray* matlabArray);
#endif
 private:

  double sigma2;
  CMatrix bias;
};

// A scaled Gaussian noise model.
class CScaleNoise : public CNoise {
 public:  
  // constructors
  CScaleNoise(){}
  CScaleNoise(const CMatrix& yin)
    {
      setTarget(yin);
      initParams();
    }
  ~CScaleNoise();
  
  double getScale(int index) const
    {
      assert(index<getNumProcesses()&&index>=0);
      return scale.getVal(index);
    }
  void setScale(double val, int index)
    {
      assert(index<getNumProcesses()&&index>=0);
      scale.setVal(val, index);
    }
  double getBias(int index) const
    {
      assert(index<getNumProcesses()&&index>=0);
      return bias.getVal(index);
    }
  void setBias(double val, int index)
    {
      assert(index<getNumProcesses()&&index>=0);
      bias.setVal(val, index);
    }
  void initStoreage();
  void initNames();
  void initVals();
  void initParams();

  ostream& display(ostream& os);
  void setParams(const CMatrix& params);
  void setParam(double val, int index);
  void getParams(CMatrix& params) const;
  double getParam(int index) const;
  void getGradParams(CMatrix& g) const;
  void getGradInputs(double& gmu, double& gvs, int i, int j) const;
  void getNuG(CMatrix& g, CMatrix& nu, int index) const;
  void updateSites(CMatrix& m, CMatrix& beta, int actIndex, const CMatrix& g, const CMatrix& nu, int index) const;
  void test(const CMatrix& muout, const CMatrix& varSigmaOut, const CMatrix& yTest) const;
  void out(CMatrix& yPred, const CMatrix& muTest, const CMatrix& varSigmaTest) const;
  void out(CMatrix& yPred, CMatrix& errorBarOut, const CMatrix& muTest, const CMatrix& varSigmaTest) const;
  void likelihoods(CMatrix& L, const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  double logLikelihood(const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  double logLikelihood() const
    {
      return logLikelihood(mu, varSigma, y);
    }
  void  getGradX(CMatrix& gX, const CMatrix& dmu, const CMatrix& cvs);
  
  inline double getMu(int i, int j) const
    {
      return mu.getVal(i, j);
    }
  inline void setMu(double val, int i, int j) 
    {
      mu.setVal(val, i, j);
    }
  inline double getVarSigma(int i, int j) const
    {
      return varSigma.getVal(i, j);
    }
  inline void setVarSigma(double val, int i, int j) 
    {
      assert(!isnan(val));
      assert(val>=0);
      varSigma.setVal(val, i, j);
    }
  inline double getTarget(int i, int j) const
    {
      return y.getVal(i, j);
    }
  void setTarget(const CMatrix& vals)
    {
      y.deepCopy(vals);
      setNumData(vals.getRows());
      setNumProcesses(vals.getCols());
    }
#ifdef _NDLMATLAB
  // Adds parameters to the mxArray.
  void addParamToMxArray(mxArray* matlabArray) const;
  // Gets the parameters from the mxArray.
  void extractParamFromMxArray(const mxArray* matlabArray);
#endif
 private:

  double sigma2;
  CMatrix bias;
  CMatrix scale;
};

// The probit noise model often used for classification.
class CProbitNoise : public CNoise {
 public:  
  // constructors
  CProbitNoise(){}
  CProbitNoise(const CMatrix& yin)
    {
      setTarget(yin);
      initParams();
    }
  ~CProbitNoise();
  
  void initStoreage();
  void initNames();
  void initVals();
  void initParams();

  ostream& display(ostream& os);
  void setParams(const CMatrix& params);
  void setParam(double val, int index);
  void getParams(CMatrix& params) const;
  double getParam(int index) const;
  void getGradParams(CMatrix& g) const;
  void getGradInputs(double& gmu, double& gvs, int i, int j) const;
  void test(const CMatrix& muout, const CMatrix& varSigmaOut, const CMatrix& yTest) const;
  void out(CMatrix& yPred, const CMatrix& muTest, const CMatrix& varSigmaTest) const;
  void out(CMatrix& yPred, CMatrix& probOut, const CMatrix& muTest, const CMatrix& varSigmaTest) const;
  void likelihoods(CMatrix& L, const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  double logLikelihood(const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  double logLikelihood() const
    {
      return logLikelihood(mu, varSigma, y);
    }
  void  getGradX(CMatrix& gX, const CMatrix& dmu, const CMatrix& cvs);
  
  inline double getMu(int i, int j) const
    {
      return mu.getVal(i, j);
    }
  inline void setMu(double val, int i, int j) 
    {
      mu.setVal(val, i, j);
    }
  inline double getVarSigma(int i, int j) const
    {
      return varSigma.getVal(i, j);
    }
  inline void setVarSigma(double val, int i, int j) 
    {
      assert(!isnan(val));
      assert(val>=0);
      varSigma.setVal(val, i, j);
    }
  inline double getTarget(int i, int j) const
    {
      return y.getVal(i, j);
    }
  virtual void setTarget(const CMatrix& vals)
    {
      y.deepCopy(vals);
      setNumData(vals.getRows());
      setNumProcesses(vals.getCols());
    }
#ifdef _NDLMATLAB
  // Adds parameters to the mxArray.
  void addParamToMxArray(mxArray* matlabArray) const;
  // Gets the parameters from the mxArray.
  void extractParamFromMxArray(const mxArray* matlabArray);
#endif

 private:
  double sigma2;
  CMatrix bias;
};

// The null category noise model for semi-supervised learning (see Lawrence and Jordan in NIPS 2004).
class CNcnmNoise : public CNoise {
 public:  
  // constructors
  CNcnmNoise(){}
  CNcnmNoise(const CMatrix& yin)
    {
      setTarget(yin);
      initParams();
    }
  ~CNcnmNoise();
  
  void initStoreage();
  void initNames();
  void initVals();
  void initParams();

  ostream& display(ostream& os);
  void setParams(const CMatrix& params);
  void setParam(double val, int index);
  void getParams(CMatrix& params) const;
  double getParam(int index) const;
  void getGradParams(CMatrix& g) const;
  void getGradInputs(double& gmu, double& gvs, int i, int j) const;
  void test(const CMatrix& muout, const CMatrix& varSigmaOut, const CMatrix& yTest) const;
  void out(CMatrix& yPred, const CMatrix& muTest, const CMatrix& varSigmaTest) const;
  void out(CMatrix& yPred, CMatrix& probOut, const CMatrix& muTest, const CMatrix& varSigmaTest) const;
  void likelihoods(CMatrix& L, const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  double logLikelihood(const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  void writeParamsToStream(ostream& out) const;
  void readParamsFromStream(istream& in);
  double logLikelihood() const
    {
      return logLikelihood(mu, varSigma, y);
    }
  void  getGradX(CMatrix& gX, const CMatrix& dmu, const CMatrix& cvs);
  
  inline double getMu(int i, int j) const
    {
      return mu.getVal(i, j);
    }
  inline void setMu(double val, int i, int j) 
    {
      mu.setVal(val, i, j);
    }
  inline double getVarSigma(int i, int j) const
    {
      return varSigma.getVal(i, j);
    }
  inline void setVarSigma(double val, int i, int j) 
    {
      if(isnan(val))
	throw ndlexceptions::Error("varSigma is being set with value NaN.");
      assert(val>=0);
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
  inline double getTarget(int i, int j) const
    {
      return y.getVal(i, j);
    }
  void setTarget(const CMatrix& vals)
    {
      y.deepCopy(vals);
      setSplitGamma(false);
      setNumData(vals.getRows());
      setNumProcesses(vals.getCols());
    }
#ifdef _NDLMATLAB
  // Adds parameters to the mxArray.
  void addParamToMxArray(mxArray* matlabArray) const;
  // Gets the parameters from the mxArray.
  void extractParamFromMxArray(const mxArray* matlabArray);
#endif

 private:
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
  COrderedNoise(){}
  COrderedNoise(const CMatrix& yin)
    {
      setTarget(yin);
      initParams();
    }
  COrderedNoise(const CMatrix& yin, int numCts) : numCats(numCts)
    {
      setTarget(yin);
      initParams();
    }
  ~COrderedNoise();
  
  int getNumCategories() const
    {
      return numCats;
    }
  void setNumCategories(int val)
    {
      assert(val>1);
      numCats=val;
    }

  void initStoreage();
  void initNames();
  void initVals();
  void initParams();

  ostream& display(ostream& os);
  void setParams(const CMatrix& params);
  void setParam(double val, int index);
  void getParams(CMatrix& params) const;
  double getParam(int index) const;
  void getGradParams(CMatrix& g) const;
  void getGradInputs(double& gmu, double& gvs, int i, int j) const;
  void test(const CMatrix& muout, const CMatrix& varSigmaOut, const CMatrix& yTest) const;
  void out(CMatrix& yPred, const CMatrix& muTest, const CMatrix& varSigmaTest) const;
  void out(CMatrix& yPred, CMatrix& probOut, const CMatrix& muTest, const CMatrix& varSigmaTest) const;
  void likelihoods(CMatrix& L, const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  double logLikelihood(const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  void writeParamsToStream(ostream& out) const;
  void readParamsFromStream(istream& in);
  double logLikelihood() const
    {
      return logLikelihood(mu, varSigma, y);
    }
  void  getGradX(CMatrix& gX, const CMatrix& dmu, const CMatrix& cvs);
  
  inline double getMu(int i, int j) const
    {
      return mu.getVal(i, j);
    }
  inline void setMu(double val, int i, int j) 
    {
      mu.setVal(val, i, j);
    }
  inline double getVarSigma(int i, int j) const
    {
      return varSigma.getVal(i, j);
    }
  inline void setVarSigma(double val, int i, int j) 
    {
      if(isnan(val))
	throw ndlexceptions::Error("varSigma is being set with value NaN.");
      assert(val>=0);
      varSigma.setVal(val, i, j);
    }
  inline double getTarget(int i, int j) const
    {
      return y.getVal(i, j);
    }
  void setTarget(const CMatrix& vals)
    {
      y.deepCopy(vals);
      setNumData(vals.getRows());
      setNumProcesses(vals.getCols());
    }
#ifdef _NDLMATLAB
  // Adds parameters to the mxArray.
  void addParamToMxArray(mxArray* matlabArray) const;
  // Gets the parameters from the mxArray.
  void extractParamFromMxArray(const mxArray* matlabArray);
#endif

 private:
  int numCats;
  CMatrix widths;
  double sigma2;
  CMatrix bias;

  mutable CMatrix gwidth;
};


void writeNoiseToStream(const CNoise& noise, ostream& out);
CNoise* readNoiseFromStream(istream& in);


#endif
