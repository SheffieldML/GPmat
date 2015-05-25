#ifndef CTRANSFORM_H
#define CTRANSFORM_H

#include <cmath>
#include "CMatrix.h"
#include "ndlutil.h"

using namespace std;

const double limVal=36;

class CTransform {

 public:
  CTransform()
    {
      transform = 0;
    }
  virtual ~CTransform()
    {
    }
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

 private:
  string type;
 protected:
  bool transform;
  const static double eps = 1e-16;
};


class CNegLogLogitTransform : public CTransform  {
 public:
  CNegLogLogitTransform();
  double atox(double a) const;
  double xtoa(double x) const;
  double gradfact(double x) const;
  
 private:
  bool transform;

};
class CLinearTransform : public CTransform  {
 public:
  CLinearTransform()
    {
      transform = 1;
      setType("linear");
      m = 1.0;
      c = 0.0;
    }
  
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
  void setM(const double val)
    {
      m = val;
    }
  double getM() const
    {
      return m;
    }
  void setC(const double val)
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
class CSigmoidTransform : public CTransform  {
 public:
  CSigmoidTransform();  
  double atox(double a) const;
  double xtoa(double x) const;
  double gradfact(double x) const;
  
 private:
  bool transform;
  
};

class CParamTransforms : CMatinterface {
  
 public:
  mxArray* toMxArray() const;
  void fromMxArray(const mxArray* transformArray);
  void addTransform(CTransform* trans, const int index)
    {
      assert(index>=0);
      transIndex.push_back(index);
      transforms.push_back(trans);
    }

  void clearTransforms()
    {
      transIndex.clear();
      transforms.clear();
    }
  inline string getTransformType(const int ind) const
    {
      assert(ind>=0);
      assert(ind<getNumTransforms());
      return transforms[ind]->getType();
    }
  inline int getTransformIndex(const int ind) const
    {
      assert(ind>=0);
      assert(ind<getNumTransforms());
      return transIndex[ind];
    }
  inline const int getNumTransforms() const
    {
      return transforms.size();
    }
  vector<CTransform*> transforms;
  vector<int> transIndex;

};





class CTransformable {

 public:

  // these are the pure virtual functions.
  virtual const int getNumParams() const=0;
  virtual double getParam(const int paramNo) const=0;
  virtual void setParam(const double val, const int paramNo)=0;
  virtual void getGradParams(CMatrix& g) const=0;

  // these are default implementations.
  virtual void getParams(CMatrix& params) const
    {
      assert(params.getRows()==1);
      assert(params.getCols()==getNumParams());
      for(int i=0; i<params.getCols(); i++)
	params.setVal(getParam(i), i);
    }
  virtual void setParams(const CMatrix& params)
    {
      assert(params.getRows()==1);
      assert(params.getCols()==getNumParams());
      for(int i=0; i<params.getCols(); i++)
	setParam(params.getVal(i), i);
    }
  
  virtual double getTransParam(const int paramNo) const
    {
      assert(paramNo>=0);
      assert(paramNo<getNumParams());
      double param = getParam(paramNo);
      vector<int>::const_iterator pos = find(transArray.transIndex.begin(), transArray.transIndex.end(), paramNo);
      if(pos == transArray.transIndex.end())
	return param;
      else
	{
	  int ind = pos - transArray.transIndex.begin();
	  return transArray.transforms[ind]->xtoa(param);
	}
    }
  virtual void getTransParams(CMatrix& transParam) const
    {
      assert(transParam.getRows()==1);
      assert(transParam.getCols()==getNumParams());
      getParams(transParam);
      double val;
      for(int i=0; i<transArray.transIndex.size(); i++)
	{
	  val=transParam.getVal(transArray.transIndex[i]);
	  transParam.setVal(transArray.transforms[i]->xtoa(val), transArray.transIndex[i]);
	}  
    }
  virtual void setTransParam(const double val, const int paramNo)
    {
      assert(paramNo>=0);
      assert(paramNo<getNumParams());
      vector<int>::iterator pos=find(transArray.transIndex.begin(), transArray.transIndex.end(), paramNo);
      if(pos==transArray.transIndex.end())
	setParam(val, paramNo);
      else
	{
	  int ind = pos - transArray.transIndex.begin();
	  setParam(transArray.transforms[ind]->atox(val), paramNo);
	}
    }
  virtual void setTransParams(const CMatrix& transParam)
    {
      assert(transParam.getRows()==1);
      assert(transParam.getCols()==getNumParams());
      CMatrix param(transParam);
      double val = 0.0;
      for(int i=0; i<transArray.transIndex.size(); i++)
	{
	  val = param.getVal(transArray.transIndex[i]);
	  param.setVal(transArray.transforms[i]->atox(val), transArray.transIndex[i]);
	}
      setParams(param);
    }
  virtual void getGradTransParams(CMatrix& g) const
    {
      assert(g.getRows()==1);
      assert(g.getCols()==getNumParams());
      getGradParams(g);
      double val;
      double param;
      for(int i=0; i<transArray.transIndex.size(); i++)
	{
	  val=g.getVal(transArray.transIndex[i]);
	  param=getParam(transArray.transIndex[i]);
	  g.setVal(val*transArray.transforms[i]->gradfact(param), transArray.transIndex[i]);
	}  
    }

  // These are non-modifiable methods.
  inline const int getNumTransforms() const
    {
      return transArray.getNumTransforms();
    }
  inline CTransform* getTransform(const int ind) const
    {
      assert(ind>=0);
      assert(ind<getNumTransforms());
      return transArray.transforms[ind];
    }
  inline string getTransformType(const int ind) const
    {
      return transArray.getTransformType(ind);
    }
  inline int getTransformIndex(const int ind) const
    {
      return transArray.getTransformIndex(ind);
    }
  inline double getTransformGradFact(const double val, const int ind) const
    {
      return transArray.transforms[ind]->gradfact(val);
    }
  void addTransform(CTransform* trans, const int index)
    {
      assert(index>=0);
      assert(index<getNumParams());
      transArray.transIndex.push_back(index);
      transArray.transforms.push_back(trans);
    }

  void clearTransforms()
    {
      transArray.transIndex.clear();
      transArray.transforms.clear();
    }

  mxArray* transformsToMxArray() const
    {
      return transArray.toMxArray();
    }
  void transformsFromMxArray(const mxArray* matlabArray) 
    {
      transArray.fromMxArray(matlabArray);
    }
 private:
  CParamTransforms transArray;

};

#endif 
