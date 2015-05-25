#ifndef CTRANSFORM_H
#define CTRANSFORM_H

#include <cmath>
#include "CMatrix.h"
#include "ndlutil.h"

using namespace std;

const double limVal=36;

// This is a base class for non-linear variable transformation.
class CTransform 
{

 public:
  CTransform() 
  {
    transform = 0;
  }
  virtual ~CTransform() {}
  CTransform(const CTransform& rhs) : type(rhs.type) {}
  virtual CTransform* clone() const=0;
  virtual double atox(double x) const=0;
  virtual double xtoa(double x) const=0;
  virtual double gradfact(double x) const=0;
  virtual void setType(const string name) 
  {
    type = name;
  }
  virtual string getType() const 
  {
    return type;
  }
  static CTransform* defaultPositive();
  static CTransform* defaultZeroOne();
  static CTransform* getNewTransformPointer(const string transformType);
 
 private:
  string type;
 protected:
  bool transform;
  #ifndef SWIG 
  // swig gives an error here
  const static double eps;
  #else 
  static double eps;
  #endif
};

// This transform transforms from real to positive numbers.
class CExpTransform : public CTransform  
{
 public:
  CExpTransform();
 CExpTransform(const CExpTransform& rhs) : CTransform(rhs), transform(rhs.transform) {}
  CTransform* clone() const {return new CExpTransform(*this);}
  double atox(double a) const;
  double xtoa(double x) const;
  double gradfact(double x) const;
  
 private:
  bool transform;

};
// This transform transforms from real to positive numbers.
class CNegLogLogitTransform : public CTransform  
{
 public:
  CNegLogLogitTransform();
 CNegLogLogitTransform(const CNegLogLogitTransform& rhs) : CTransform(rhs), transform(rhs.transform) {}
  CTransform* clone() const {return new CNegLogLogitTransform(*this);}
  double atox(double a) const;
  double xtoa(double x) const;
  double gradfact(double x) const;
  
 private:
  bool transform;

};
class CLinearTransform : public CTransform  
{
 public:
  CLinearTransform()
  {
    transform = 1;
    setType("linear");
    m = 1.0;
    c = 0.0;
  }
 CLinearTransform(const CLinearTransform& rhs) : CTransform(rhs), transform(rhs.transform), m(rhs.m), c(rhs.c) {}
  CTransform* clone() const {return new CLinearTransform(*this);}
  
  double atox(double a) const
  {
    return (a-c)/m;
  }
  double xtoa(double x) const 
  {
    return m*x+c;
  }
  double gradfact(double x) const
  {
    return 1/m;
  }
  void setM(double val)
  {
    m = val;
  }
  double getM() const
  {
    return m;
  }
  void setC(double val)
  {
    c = val;
  }
  double getC() const
  {
    return c;
  }
 private:
  bool transform;
  double m;
  double c;

};

// This transformation goes from real numbers to the range [0->1]
class CSigmoidTransform : public CTransform  
{
 public:
  CSigmoidTransform();  
 CSigmoidTransform(const CSigmoidTransform& rhs) : CTransform(rhs), transform(rhs.transform) {}
  CTransform* clone() const {return new CSigmoidTransform(*this);}
  double atox(double a) const;
  double xtoa(double x) const;
  double gradfact(double x) const;
  
 private:
  bool transform;
  
};

// A class for storing the parameter transformations.
class CParamTransforms : public CMatInterface, public CStreamInterface
{
  
 public:
  string getType() const
  { 
    return "transforms";
  }
  string getBaseType() const
  {
    return "transforms";
  }
  bool equals(CParamTransforms transforms) const;
  void display(ostream& out) const;
  void writeParamsToStream(ostream& out) const;
  void readParamsFromStream(istream& in);
#ifdef _NDLMATLAB
  mxArray* toMxArray() const;
  void fromMxArray(const mxArray* transformArray);
#endif
  void addTransform(const CTransform* trans, unsigned int index) 
  {
    
    transIndex.push_back(index);
    transforms.push_back(trans->clone());
  }

  void clearTransforms() 
  {
    transIndex.clear();
    transforms.clear();
  }
  inline string getTransformType(unsigned int ind) const 
  {
    
    BOUNDCHECK(ind<getNumTransforms());
    return transforms[ind]->getType();
  }
  inline unsigned int getTransformIndex(unsigned int ind) const 
  {
    
    BOUNDCHECK(ind<getNumTransforms());
    return transIndex[ind];
  }
  inline unsigned int getNumTransforms() const 
  {
    return transforms.size();
  }
	

  vector<CTransform*> transforms;
  vector<unsigned int> transIndex;
};




// This is an abstract base class for making the parameters of a class transformable.
class CTransformable 
{

 public:

  // these are the pure virtual functions.
  virtual ~CTransformable(){}
  virtual unsigned int getNumParams() const=0;
  virtual double getParam(unsigned int paramNo) const=0;
  virtual void setParam(double val, unsigned int paramNo)=0;
  virtual void getGradParams(CMatrix& g) const=0;

  // these are default implementations.
  virtual void getParams(CMatrix& params) const 
  {
    if(params.getRows()!=1)
      throw ndlexceptions::RuntimeError("getParams(): Dimension match check failed, numbers of rows should be 1, currently it is " + ndlstrutil::itoa(params.getRows()));
    if(params.getCols()!=getNumParams())
      throw ndlexceptions::RuntimeError("getParams(): Dimension match check failed, numbers of columns should be " + ndlstrutil::itoa(getNumParams()) + ", currently it is " + ndlstrutil::itoa(params.getRows()));
    for(unsigned int i=0; i<params.getCols(); i++)
      params.setVal(getParam(i), i);
  }
  virtual void setParams(const CMatrix& params) 
  {
    if(params.getRows()!=1)
      throw ndlexceptions::RuntimeError("setParams(): Dimension match check failed, numbers of rows should be 1, currently it is " + ndlstrutil::itoa(params.getRows()));
    if(params.getCols()!=getNumParams())
      throw ndlexceptions::RuntimeError("setParams(): Dimension match check failed, numbers of columns should be " + ndlstrutil::itoa(getNumParams()) + ", currently it is " + ndlstrutil::itoa(params.getRows()));
    for(unsigned int i=0; i<params.getCols(); i++)
      setParam(params.getVal(i), i);
  }
  
  virtual double getTransParam(unsigned int paramNo) const 
  {
    BOUNDCHECK(paramNo<getNumParams());
    double param = getParam(paramNo);
    vector<unsigned int>::const_iterator pos = find(transArray.transIndex.begin(), 
					   transArray.transIndex.end(), 
					   paramNo);
    if(pos == transArray.transIndex.end())
      return param;
    else 
    {
      unsigned int ind = pos - transArray.transIndex.begin();
      return transArray.transforms[ind]->xtoa(param);
    }
  }
  virtual void getTransParams(CMatrix& transParam) const 
  {
    if(transParam.getRows()!=1)
      throw ndlexceptions::RuntimeError("getTransParams(): Dimension match check failed, numbers of rows should be 1, currently it is " + ndlstrutil::itoa(transParam.getRows()));
    if(transParam.getCols()!=getNumParams())
      throw ndlexceptions::RuntimeError("getTransParams(): Dimension match check failed, numbers of columns should be " + ndlstrutil::itoa(getNumParams()) + ", currently it is " + ndlstrutil::itoa(transParam.getRows()));
    getParams(transParam);
    double val;
    for(unsigned int i=0; i<transArray.transIndex.size(); i++) {
      val=transParam.getVal(transArray.transIndex[i]);
      transParam.setVal(transArray.transforms[i]->xtoa(val), transArray.transIndex[i]);
    }  
  }
  virtual void setTransParam(double val, unsigned int paramNo) 
  {
    BOUNDCHECK(paramNo<getNumParams());
    // this casting is required under solaris for some reason
    vector<unsigned int>::iterator pos=find(transArray.transIndex.begin(), 
					    transArray.transIndex.end(), 
					    paramNo);
    if(pos==transArray.transIndex.end())
      setParam(val, paramNo);
    else {
      unsigned int ind = pos - transArray.transIndex.begin();
      setParam(transArray.transforms[ind]->atox(val), paramNo);
    }
  }
  virtual void setTransParams(const CMatrix& transParam) 
  {
    if(transParam.getRows()!=1)
      throw ndlexceptions::RuntimeError("setTransParams(): Dimension match check failed, numbers of rows should be 1, currently it is " + ndlstrutil::itoa(transParam.getRows()));
    if(transParam.getCols()!=getNumParams())
      throw ndlexceptions::RuntimeError("setTransParams(): Dimension match check failed, numbers of columns should be " + ndlstrutil::itoa(getNumParams()) + ", currently it is " + ndlstrutil::itoa(transParam.getRows()));
    CMatrix param(transParam);
    double val = 0.0;
    for(unsigned int i=0; i<transArray.transIndex.size(); i++) 
    {
      val = param.getVal(transArray.transIndex[i]);
      param.setVal(transArray.transforms[i]->atox(val), transArray.transIndex[i]);
    }
    setParams(param);
  }
  virtual void getGradTransParams(CMatrix& g) const 
  {
    if(g.getRows()!=1)
      throw ndlexceptions::RuntimeError("getGradTransParams(): Dimension match check failed, numbers of rows should be 1, currently it is " + ndlstrutil::itoa(g.getRows()));
    if(g.getCols()!=getNumParams())
      throw ndlexceptions::RuntimeError("getGradTransParams(): Dimension match check failed, numbers of columns should be " + ndlstrutil::itoa(getNumParams()) + ", currently it is " + ndlstrutil::itoa(g.getRows()));
    getGradParams(g);
    double val;
    double param;
    for(size_t i=0; i<transArray.transIndex.size(); i++) 
    {
      val=g.getVal(transArray.transIndex[i]);
      param=getParam(transArray.transIndex[i]);
      g.setVal(val*transArray.transforms[i]->gradfact(param), transArray.transIndex[i]);
    }  
  }
  
  // These are non-modifiable methods.
  inline unsigned int getNumTransforms() const 
  {
    return transArray.getNumTransforms();
  }
  inline const CTransform* getTransform(unsigned int ind) const 
  {
    
    BOUNDCHECK(ind<getNumTransforms());
    return transArray.transforms[ind];
  }
  inline string getTransformType(unsigned int ind) const 
  {
    return transArray.getTransformType(ind);
  }
  inline unsigned int getTransformIndex(unsigned int ind) const 
  {
    return transArray.getTransformIndex(ind);
  }
  inline double getTransformGradFact(double val, unsigned int ind) const 
  {
    return transArray.transforms[ind]->gradfact(val);
  }
  void addTransform(const CTransform* trans, unsigned int index) 
  {
    
    BOUNDCHECK(index<getNumParams());
    transArray.transIndex.push_back(index);
    transArray.transforms.push_back(trans->clone());
  }
  
  void clearTransforms() 
  {
    for(size_t i = 0; i<transArray.transforms.size(); i++)
      delete transArray.transforms[i];
    transArray.transIndex.clear();
    transArray.transforms.clear();
  }
  
#ifdef _NDLMATLAB
  mxArray* transformsToMxArray() const 
  {
    return transArray.toMxArray();
  }
  void transformsFromMxArray(const mxArray* matlabArray)  
  {
    transArray.fromMxArray(matlabArray);
  }
#endif
 private:
  CParamTransforms transArray;

};


#endif 

