#ifndef CNOISE_H
#define CNOISE_H
#include <iostream>
#include <cmath>
#include <vector>
#include "CMatrix.h"
#include "CTransform.h"
#include "CDist.h"
#include "CKern.h"
using namespace std;

class CNoise {

 public:
  // constructors
  CNoise(){}
  CNoise(const CMatrix& y){}
  CNoise(mxArray* noise){}
  virtual ~CNoise(){}
  CNoise(const CNoise& noise) {}
  
  // pure virtual functions
  virtual void setInitParam()=0;
  virtual void setInitParam(const CMatrix& y)=0;
  virtual ostream& display(ostream& os)=0;
  virtual void setParams(const CMatrix& X)=0;
  virtual void getParams(CMatrix& X) const=0;
  virtual void getGradParams(CMatrix& g, const CMatrix& mu, const CMatrix& varsigma, const CMatrix& y) const=0;
  virtual void getGradInputs(CMatrix& dlnZ_dmu, CMatrix& dlnZ_dvs, const CMatrix& mu, const CMatrix& varsigma, const CMatrix& y) const=0;
  virtual void getNuG(CMatrix& g, CMatrix& nu, const CMatrix& mu, const CMatrix& varsigma, const CMatrix& y, const int index) const=0;
  virtual void updateSites(CMatrix& m, CMatrix& beta, const int actIndex, const CMatrix& g, const CMatrix& nu, const CMatrix& mu, const CMatrix& varsigma, const CMatrix& y, const int index) const=0;
  virtual CMatrix out(const CMatrix& mu, const CMatrix& varsigma) const=0;
  virtual CMatrix likelihood(const CMatrix& mu, const CMatrix& varsigma, const CMatrix& y) const=0;
  virtual double logLikelihood(const CMatrix& mu, const CMatrix& varsigma, const CMatrix& y) const=0;
  
  // non virtual functions
  void getGradX(CMatrix& gX, const CMatrix& mu, const CMatrix& varsigma, const CMatrix& dmu, const CMatrix& cvs, const CMatrix& y);
  void getGradTransParams(CMatrix& g, const CMatrix& mu, const CMatrix& varsigma, const CMatrix& y);
  void getTransParams(CMatrix& params) const;
  void setTransParams(const CMatrix& params);
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
  inline int getNParams() const
    {
      return nParams;
    }
  inline int getNProcess() const
    {
      return nProcess;
    }
  void addTransform(CTransform* trans, const int index)
    {
      assert(index<nParams);
      transIndex.push_back(index);
      transforms.push_back(trans);
    }
 protected:
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
  virtual inline void setNProcess(const int num)
    {
      nProcess = num;
    }
  inline void setNParams(const int num) 
    {
      nParams = num;
    }
  
  vector<CDist*> priors;
  bool spherical;
  bool logConcave;
  bool missing;
  int nParams;
  int nProcess;
  vector<CTransform*> transforms;
  vector<int> transIndex;
};

class CGaussianNoise : public CNoise {
 public:
  // constructors
  CGaussianNoise(){}
  CGaussianNoise(const CMatrix& y)
    {
      bias.resize(1, y.getCols());
      setInitParam(y);
    }
  CGaussianNoise(mxArray* noise){}
  ~CGaussianNoise();
  CGaussianNoise(const CGaussianNoise& noise){}
  
  void setInitParam();
  void setInitParam(const int);
  void setInitParam(const CMatrix& y);
  ostream& display(ostream& os);
  void setParams(const CMatrix& params);
  void getParams(CMatrix& params) const;
  void getGradParams(CMatrix& g, const CMatrix& mu, const CMatrix& varsigma, const CMatrix& y) const;
  void getGradInputs(CMatrix& dlnZ_dmu, CMatrix& dlnZ_dvs, const CMatrix& mu, const CMatrix& varsigma, const CMatrix& y) const;
  void getNuG(CMatrix& g, CMatrix& nu, const CMatrix& mu, const CMatrix& varsigma, const CMatrix& y, const int index) const;
  void updateSites(CMatrix& m, CMatrix& beta, const int actIndex, const CMatrix& g, const CMatrix& nu, const CMatrix& mu, const CMatrix& varsigma, const CMatrix& y, const int index) const;
  CMatrix out(const CMatrix& mu, const CMatrix& varsigma) const;
  CMatrix likelihood(const CMatrix& mu, const CMatrix& varsigma, const CMatrix& y) const;
  double logLikelihood(const CMatrix& mu, const CMatrix& varsigma, const CMatrix& y)const;
  void  getGradX(CMatrix& gX, const CMatrix& mu, const CMatrix& varsigma, const CMatrix& dmu, const CMatrix& cvs, const CMatrix& y);
  
 private:
  double sigma2;
  CMatrix bias;
};

#endif
