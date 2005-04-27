#ifndef CNOISE_H
#define CNOISE_H
#include <iostream>
#include <cmath>
#include <climits>
#include "ndlutil.h"
#include <vector>
#include "CMatrix.h"
#include "CTransform.h"
#include "COptimisable.h"
#include "CDist.h"
#include "CKern.h"
using namespace std;

class CNoise : public CTransformable, public COptimisable, public CMatinterface {

 public:
  // constructors
  CNoise() {}
  virtual ~CNoise(){}
  
  // pure virtual functions
  virtual void setInitParam()=0;
  virtual ostream& display(ostream& os)=0;
  virtual void setParams(const CMatrix& X)=0;
  virtual void getParams(CMatrix& X) const=0;
  virtual void getGradParams(CMatrix& g) const=0;
  virtual void getGradInputs(double& dlnZ_dmu, double& dlnZ_dvs, const int i, const int j) const=0;
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
  virtual void getNuG(CMatrix& g, CMatrix& nu, const int index) const;
  virtual void updateSites(CMatrix& m, CMatrix& beta, const int actIndex, const CMatrix& g, const CMatrix& nu, const int index) const;
  virtual void out(CMatrix& yTest, const CMatrix& muTest, const CMatrix& varSigmaTest) const=0;
  virtual void likelihood(CMatrix& L, const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const=0;
  virtual double logLikelihood(const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const=0;
  virtual double logLikelihood() const=0;
  virtual mxArray* toMxArray() const;
  virtual void fromMxArray(const mxArray* matlabArray);

  // Adds parameters to the mxArray.
  virtual void addParamToMxArray(mxArray* matlabArray) const;
  // Gets the parameters from the mxArray.
  virtual void extractParamFromMxArray(const mxArray* matlabArray);

  // non virtual functions
  inline void setVerbosity(const int val)
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
  inline const int getNumParams() const
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
  virtual string getParamName(const int index) const
    {
      assert(index>=0);
      assert(index<paramNames.size());
      return paramNames[index];
    }
  void setParamName(const string paramName, const int index)
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

  virtual void setMu(const double val, const int i, const int j)=0;
  virtual double getMu(const int i, const int j) const=0;
  virtual void setVarSigma(const double val, const int i, const int j)=0;
  virtual double getVarSigma(const int i, const int j) const=0;
  virtual double getTarget(const int i, const int j) const=0;
  virtual void setTarget(const CMatrix& vals)=0;

  virtual void setMus(const double val)
    {
      for(int i=0; i<getNumData(); i++)
	for(int j=0; j<getNumProcesses(); j++)
	  setMu(val, i, j);
    }
  virtual void setVarSigmas(const double val)
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

  bool equals(const CNoise& noise, const double tol=ndlutil::MATCHTOL) const;
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
  inline void setNumParams(const int num) 
    {
      nParams = num;
    }
  inline void setNumData(const int num) 
    {
      nData = num;
    }
  inline void setNumProcesses(const int num) 
    {
      nProcesses = num;
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

class CGaussianNoise : public CNoise {
 public:  
  // constructors
  CGaussianNoise(){}
  CGaussianNoise(const CMatrix& yin)
    {
      setTarget(yin);
    }
  ~CGaussianNoise();
  
  void setInitParam();
  ostream& display(ostream& os);
  void setParams(const CMatrix& params);
  void setParam(const double val, const int index);
  void getParams(CMatrix& params) const;
  double getParam(const int index) const;
  void getGradParams(CMatrix& g) const;
  void getGradInputs(double& gmu, double& gvs, const int i, const int j) const;
  void getNuG(CMatrix& g, CMatrix& nu, const int index) const;
  void updateSites(CMatrix& m, CMatrix& beta, const int actIndex, const CMatrix& g, const CMatrix& nu, const int index) const;
  void out(CMatrix& yTest, const CMatrix& muTest, const CMatrix& varSigmaTest) const;
  void likelihood(CMatrix& L, const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  double logLikelihood(const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  double logLikelihood() const
    {
      return logLikelihood(mu, varSigma, y);
    }
  void  getGradX(CMatrix& gX, const CMatrix& dmu, const CMatrix& cvs);
  
  inline double getMu(const int i, const int j) const
    {
      return mu.getVal(i, j);
    }
  inline void setMu(const double val, const int i, const int j) 
    {
      mu.setVal(val, i, j);
    }
  inline double getVarSigma(const int i, const int j) const
    {
      return varSigma.getVal(i, j);
    }
  inline void setVarSigma(const double val, const int i, const int j) 
    {
      assert(!isnan(val));
      assert(val>=0);
      varSigma.setVal(val, i, j);
    }
  inline double getTarget(const int i, const int j) const
    {
      return y.getVal(i, j);
    }
  virtual void setTarget(const CMatrix& vals)
    {
      y.deepCopy(vals);
      setNumData(vals.getRows());
      setNumProcesses(vals.getCols());
      setInitParam();
    }
  // Adds parameters to the mxArray.
  void addParamToMxArray(mxArray* matlabArray) const;
  // Gets the parameters from the mxArray.
  void extractParamFromMxArray(const mxArray* matlabArray);

 private:

  double sigma2;
  CMatrix bias;
};


class CProbitNoise : public CNoise {
 public:  
  // constructors
  CProbitNoise(){}
  CProbitNoise(const CMatrix& yin)
    {
      setTarget(yin);
    }
  ~CProbitNoise();
  
  void setInitParam();
  ostream& display(ostream& os);
  void setParams(const CMatrix& params);
  void setParam(const double val, const int index);
  void getParams(CMatrix& params) const;
  double getParam(const int index) const;
  void getGradParams(CMatrix& g) const;
  void getGradInputs(double& gmu, double& gvs, const int i, const int j) const;
  void out(CMatrix& yTest, const CMatrix& muTest, const CMatrix& varSigmaTest) const;
  void likelihood(CMatrix& L, const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  double logLikelihood(const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  double logLikelihood() const
    {
      return logLikelihood(mu, varSigma, y);
    }
  void  getGradX(CMatrix& gX, const CMatrix& dmu, const CMatrix& cvs);
  
  inline double getMu(const int i, const int j) const
    {
      return mu.getVal(i, j);
    }
  inline void setMu(const double val, const int i, const int j) 
    {
      mu.setVal(val, i, j);
    }
  inline double getVarSigma(const int i, const int j) const
    {
      return varSigma.getVal(i, j);
    }
  inline void setVarSigma(const double val, const int i, const int j) 
    {
      assert(!isnan(val));
      assert(val>=0);
      varSigma.setVal(val, i, j);
    }
  inline double getTarget(const int i, const int j) const
    {
      return y.getVal(i, j);
    }
  virtual void setTarget(const CMatrix& vals)
    {
      y.deepCopy(vals);
      setNumData(vals.getRows());
      setNumProcesses(vals.getCols());
      setInitParam();
    }
  // Adds parameters to the mxArray.
  void addParamToMxArray(mxArray* matlabArray) const;
  // Gets the parameters from the mxArray.
  void extractParamFromMxArray(const mxArray* matlabArray);
  

 private:
  double sigma2;
  CMatrix bias;
};

class CNcnmNoise : public CNoise {
 public:  
  // constructors
  CNcnmNoise(){}
  CNcnmNoise(const CMatrix& yin)
    {
      setTarget(yin);
    }
  ~CNcnmNoise();
  
  void setInitParam();
  ostream& display(ostream& os);
  void setParams(const CMatrix& params);
  void setParam(const double val, const int index);
  void getParams(CMatrix& params) const;
  double getParam(const int index) const;
  void getGradParams(CMatrix& g) const;
  void getGradInputs(double& gmu, double& gvs, const int i, const int j) const;
  void out(CMatrix& yTest, const CMatrix& muTest, const CMatrix& varSigmaTest) const;
  void likelihood(CMatrix& L, const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  double logLikelihood(const CMatrix& muTest, const CMatrix& varSigmaTest, const CMatrix& yTest) const;
  double logLikelihood() const
    {
      return logLikelihood(mu, varSigma, y);
    }
  void  getGradX(CMatrix& gX, const CMatrix& dmu, const CMatrix& cvs);
  
  inline double getMu(const int i, const int j) const
    {
      return mu.getVal(i, j);
    }
  inline void setMu(const double val, const int i, const int j) 
    {
      mu.setVal(val, i, j);
    }
  inline double getVarSigma(const int i, const int j) const
    {
      return varSigma.getVal(i, j);
    }
  inline void setVarSigma(const double val, const int i, const int j) 
    {
      if(isnan(val))
	throw("varSigma is being set with value NaN.");
      assert(val>=0);
      varSigma.setVal(val, i, j);
    }
  inline double getTarget(const int i, const int j) const
    {
      return y.getVal(i, j);
    }
  virtual void setTarget(const CMatrix& vals)
    {
      y.deepCopy(vals);
      setNumData(vals.getRows());
      setNumProcesses(vals.getCols());
      setInitParam();
    }
  // Adds parameters to the mxArray.
  void addParamToMxArray(mxArray* matlabArray) const;
  // Gets the parameters from the mxArray.
  void extractParamFromMxArray(const mxArray* matlabArray);
  

 private:
  double gamman;
  double gammap;
  double width;
  double sigma2;
  CMatrix bias;
};

#endif
