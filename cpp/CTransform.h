#ifndef CTRANSFORM_H
#define CTRANSFORM_H

#include <cmath>
#include "CMatrix.h"

using namespace std;
class CTransform {

 public:
  CTransform()
    {
      transform = 0;
    }
  virtual ~CTransform()
    {
    }
  virtual double atox(double x)=0;
  virtual double xtoa(double x)=0;
  virtual double gradfact(double x)=0;
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
  const static double limVal = 36;
};


class CNegLogLogitTransform : public CTransform  {
 public:
  CNegLogLogitTransform()
    {
      transform = 1;
      setType("negLogLogit");
    }
  
  double atox(double a)
    {
      double x;
      if(a<-limVal)
	x = exp(-limVal);
      else if(a<limVal)
	x=log(1+exp(a));
      return x;
    }
  double xtoa(double x)
    {
      double a=x;
      assert(a>-limVal);
      if(a<limVal)
	a=log(exp(x)-1);
      return a;
    }
  double gradfact(double x)
    {
      double g;
      assert(x>-limVal);
      if(x<limVal)
	g=(exp(x)-1)/exp(x);
      else
	g=1.0;
      return g;
    }

  
 private:
  bool transform;
  const static double limVal = 36;

};

class CTransformable {

 public:
  virtual const int getNumParams() const=0;
  virtual double getParam(const int paramNo) const=0;
  virtual void setParam(const double val, const int paramNo)=0;
  virtual void getParams(CMatrix& params) const=0;
  virtual void setParams(const CMatrix& params)=0;
  virtual void getGradParams(CMatrix& g) const=0;
  
  virtual double getTransParam(const int paramNo) const
    {
      assert(paramNo>=0);
      assert(paramNo<getNumParams());
      double param = getParam(paramNo);
      vector<int>::const_iterator pos = find(transIndex.begin(), transIndex.end(), paramNo);
      if(pos == transIndex.end())
	return param;
      else
	{
	  int ind = pos - transIndex.begin();
	  return transforms[ind]->xtoa(param);
	}
    }
  virtual void getTransParams(CMatrix& transParam) const
    {
      assert(transParam.getRows()==1);
      assert(transParam.getCols()==getNumParams());
      getParams(transParam);
      double val;
      for(int i=0; i<transIndex.size(); i++)
	{
	  val=transParam.getVals(transIndex[i]);
	  transParam.setVals(transforms[i]->xtoa(val), transIndex[i]);
	}  
    }
  virtual void setTransParam(const double val, const int paramNo)
    {
      assert(paramNo>=0);
      assert(paramNo<getNumParams());
      vector<int>::iterator pos=find(transIndex.begin(), transIndex.end(), paramNo);
      if(pos==transIndex.end())
	setParam(val, paramNo);
      else
	{
	  int ind = pos - transIndex.begin();
	  setParam(transforms[ind]->atox(val), paramNo);
	}
    }
  virtual void setTransParams(const CMatrix& transParam)
    {
      assert(transParam.getRows()==1);
      assert(transParam.getCols()==getNumParams());
      CMatrix param(transParam);
      double val = 0.0;
      for(int i=0; i<transIndex.size(); i++)
	{
	  val = param.getVals(transIndex[i]);
	  param.setVals(transforms[i]->atox(val), transIndex[i]);
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
      for(int i=0; i<transIndex.size(); i++)
	{
	  val=g.getVals(transIndex[i]);
	  param=getParam(transIndex[i]);
	  g.setVals(val*transforms[i]->gradfact(param), transIndex[i]);
	}  
    }
  inline const int getNumTransforms() const
    {
      return transforms.size();
    }
  inline CTransform* getTransform(const int ind) const
    {
      assert(ind>=0);
      assert(ind<getNumTransforms());
      return transforms[ind];
    }
  inline string getTransformType(const int ind) const
    {
      return transforms[ind]->getType();
    }
  inline int getTransformIndex(const int ind) const
    {
      assert(ind>=0);
      assert(ind<getNumTransforms());
      return transIndex[ind];
    }
  inline double getTransformGradFact(const double val, const int ind) const
    {
      return transforms[ind]->gradfact(val);
    }
  void addTransform(CTransform* trans, const int index)
    {
      assert(index>=0);
      assert(index<getNumParams());
      transIndex.push_back(index);
      transforms.push_back(trans);
    }

 private:
  vector<CTransform*> transforms;
  vector<int> transIndex;

};
#endif 
