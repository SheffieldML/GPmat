#ifndef CMLTOOLS_H
#define CMLTOOLS_H

#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "ndlassert.h"
#include "ndlexceptions.h"
#include "ndlstrutil.h"
#include "COptimisable.h"
#include "CNdlInterfaces.h"
#include "CMatrix.h"
#include "CKern.h"
#include "CNoise.h"
#include "CTransform.h"
#include "CMatrix.h"
#include "CDist.h"
#include "CDataModel.h"
#include "ndlstrutil.h"
#include "ndlutil.h"

using namespace std;

const string MLTOOLSVERSION="0.1";



// Multi-layer peceptron model.
class CLinearMapping: public CMapModel, public COptimisable, public CStreamInterface, public CMatInterface
{
 public:
  CLinearMapping();
  CLinearMapping(CMatrix* pXin, CMatrix* pyOut, int verbos=2);

  //CLinearMapping(const CMatrix* X);
  ~CLinearMapping(){};
  //CLinearMapping(const CLinearMapping& model);
  CLinearMapping * clone() const
  {
    return new CLinearMapping (*this);
  }


#ifdef _NDLMATLAB
  // Constructor using file containing linearMapping.
  CLinearMapping(CMatrix* inData, 
		 CMatrix* targetData, 
		 const string linearMappingInfoFile, 
		 const string linearMappingInfoVariable, 
		 int verbos=2);
  mxArray* toMxArray() const;
  void fromMxArray(const mxArray* matlabArray);
  
#endif
  
  void initStoreage();
  void initVals();
  //void setInitParam();
  //void setParam(double val, unsigned int paramNum);
  //double getParam(int paramNum) const;
  //double getGradParam(int index, const CMatrix& X) const;
  
  void out(CMatrix& yPred, const CMatrix& inData) const;
  //void out(CMatrix& yPred, CMatrix& probPred, const CMatrix& inData) const;
  double outGradParams(CMatrix& g, const CMatrix& inData, unsigned int i, unsigned int j) const;
  double outGradX(CMatrix& g, const CMatrix& inData, unsigned int i, unsigned int j) const;
  
  void setWeights(const CMatrix& W, unsigned int layer);
  void setBias(const CMatrix& b, unsigned int layer);
  
  
  double logLikelihood() const;
  double logLikelihoodGradient(CMatrix& g) const;
  void updateG() const;
  double pointLogLikelihood(const CMatrix& yOut, const CMatrix& Xin) const;
  void display(ostream& os) const;
  
  void writeParamsToStream(ostream& os) const;
  void readParamsFromStream(istream& is);
  

  void optimise(unsigned int iters=1000);
  // specify tests for equality between two models.
  virtual bool equals(const CLinearMapping& model, double tol=ndlutil::MATCHTOL) const;
  
  virtual unsigned int getOptNumParams() const;
  virtual void getOptParams(CMatrix& param) const;
  virtual void setOptParams(const CMatrix& param);
  
  double computeObjectiveGradParams(CMatrix& g) const 
  {
    double L = logLikelihoodGradient(g);
    g.negate();
    return -L;
  }
  double computeObjectiveVal() const 
  {
    return -logLikelihood();
  }
  
 private:
  CMatrix W;
  CMatrix b;
  mutable CMatrix outActive;
  double variance;
  CMatrix* pX;
  CMatrix* py;

  void _init(); // initialise the structure.
};


// Multi-layer peceptron model.
class CMlpMapping: public CMapModel, public COptimisable, public CStreamInterface, public CMatInterface
{
 public:
  CMlpMapping();
  CMlpMapping(CMatrix* pXin, CMatrix* pyOut, unsigned int nhidden, int verbos=2);

  //CMlpMapping(const CMatrix* X);
  ~CMlpMapping(){};
  //CMlpMapping(const CMlpMapping& model);
  CMlpMapping * clone() const
  {
    return new CMlpMapping (*this);
  }


#ifdef _NDLMATLAB
  // Constructor using file containing mlpMapping.
  CMlpMapping(CMatrix* inData, 
       CMatrix* targetData, 
       const string mlpMappingInfoFile, 
       const string mlpMappingInfoVariable, 
       int verbos=2);
  mxArray* toMxArray() const;
  void fromMxArray(const mxArray* matlabArray);

#endif

  void initStoreage();
  void initVals();
  //void setInitParam();
  //void setParam(double val, unsigned int paramNum);
  //double getParam(int paramNum) const;
  //double getGradParam(int index, const CMatrix& X) const;

  void out(CMatrix& yPred, const CMatrix& inData) const;
  //void out(CMatrix& yPred, CMatrix& probPred, const CMatrix& inData) const;
  double outGradParams(CMatrix& g, const CMatrix& inData, unsigned int i, unsigned int j) const;
  double outGradX(CMatrix& g, const CMatrix& inData, unsigned int i, unsigned int j) const;

  void setWeights(const CMatrix& W, unsigned int layer);
  void setBias(const CMatrix& b, unsigned int layer);


  double logLikelihood() const;
  double logLikelihoodGradient(CMatrix& g) const;
  void updateG() const;
  double pointLogLikelihood(const CMatrix& yOut, const CMatrix& Xin) const;
  void display(ostream& os) const;

  void writeParamsToStream(ostream& os) const;
  void readParamsFromStream(istream& is);

  unsigned int getHiddenDim() const
  {
    return hiddenDim;
  }

  void optimise(unsigned int iters=1000);
  // specify tests for equality between two models.
  virtual bool equals(const CMlpMapping& model, double tol=ndlutil::MATCHTOL) const;

  virtual unsigned int getOptNumParams() const;
  virtual void getOptParams(CMatrix& param) const;
  virtual void setOptParams(const CMatrix& param);

  double computeObjectiveGradParams(CMatrix& g) const 
  {
    double L = logLikelihoodGradient(g);
    g.negate();
    return -L;
  }
  double computeObjectiveVal() const 
  {
    return -logLikelihood();
  }
  
 private:
  CMatrix W1;
  CMatrix b1;
  CMatrix W2;
  CMatrix b2;
  mutable CMatrix hiddenActive;
  mutable CMatrix tanhActive;
  mutable CMatrix outActive;
  double variance;
  unsigned int hiddenDim;
  CMatrix* pX;
  CMatrix* py;

  void _init(); // initialise the structure.
};

//void writeMapModelToStream(const CMapModel& model, ostream& out);
//CMapModel* readMapModelFromStream(istream& in);


#endif
