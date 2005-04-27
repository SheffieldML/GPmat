#ifndef CDIST_H
#define CDIST_H
#include <cmath>
#include <vector>
#include "CMatrix.h"
#include "CTransform.h"
#include "ndlutil.h"
#include "mex.h"


class CDist : public CTransformable {
  
 public:
  CDist(){}
  ~CDist(){}
  const int getNumParams() const
    {
      return nParams;
    }
  void setNumParams(const int num)
    {
      nParams = num;
    }
  virtual double getParam(const int paramNo) const=0;
  virtual void setParam(const double val, const int paramNo)=0;
  virtual void getGradParams(CMatrix& g) const
    {
      // This is a dummy function
      cerr << "getGradParams should not be used in CDist" << endl;
      exit(1);
    }
  //CDist(CDist& dist);
  virtual double getGradInput(double x) const=0;
  void setInitParam();

  virtual double logProb(double val) const=0;
  void setParamName(const string name, const int index)
    {
      assert(index>=0);
      assert(index<nParams);
      if(paramNames.size() == index)
	paramNames.push_back(name);
      else 
	{
	  if(paramNames.size()<index)
	    paramNames.resize(index+1, "no name");
	  paramNames[index] = name;
	}
    }
  virtual string getParamName(const int index) const
    {
      assert(index>=0);
      assert(index<paramNames.size());
      return paramNames[index];
    }

  void setType(string name)
    {
      type = name;
    }
  string getType() const
    {
      return type;
    }
  void setName(string name)
    {
      distName = name;
    }
  string getName() const
    {
      return distName;
    }
  
 private:
  int nParams;
  string type;
  string distName;
  vector<string> paramNames;
};


class CGaussianDist : public CDist {

 public:
  CGaussianDist();
  CGaussianDist(const CGaussianDist&);
  ~CGaussianDist();
  CGaussianDist* clone() const
    {
      return new CGaussianDist(*this);
    }
  double getParam(const int paramNo) const;
  void setParam(const double val, const int paramNo);
  double getGradInput(double x) const;
  void setInitParam();
  double logProb(double val) const;

 private:
  double precision;
};


class CGammaDist : public CDist {

 public:
  CGammaDist();
  CGammaDist(const CGammaDist&);
  ~CGammaDist();
  CGammaDist* clone() const
    {
      return new CGammaDist(*this);
    }
  double getParam(const int paramNo) const;
  void setParam(const double val, const int paramNo);
  double getGradInput(double x) const;
  void setInitParam();
  double logProb(double val) const;

 private:
  double a;
  double b;
};

class CParamPriors : CMatinterface {
  
 public:
  mxArray* toMxArray() const;
  void fromMxArray(const mxArray* distArray);
  void addDist(CDist* dist, const int index)
    {
      assert(index>=0);
      distIndex.push_back(index);
      dists.push_back(dist);
    }

  void clearDists()
    {
      distIndex.clear();
      dists.clear();
    }
  inline string getDistType(const int ind) const
    {
      assert(ind>=0);
      assert(ind<getNumDists());
      return dists[ind]->getType();
    }
  inline int getDistIndex(const int ind) const
    {
      assert(ind>=0);
      assert(ind<getNumDists());
      return distIndex[ind];
    }
  inline const int getNumDists() const
    {
      return dists.size();
    }
  vector<CDist*> dists;
  vector<int> distIndex;

};


class CRegularisable {

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
  
  virtual void addPriorGrad(CMatrix& g) const
    {
      assert(g.getRows()==1);
      assert(g.getCols()==getNumParams());
      double val;
      double param;
      for(int i=0; i<distArray.distIndex.size(); i++)
	{
	  val=g.getVal(distArray.distIndex[i]);
	  param=getParam(distArray.distIndex[i]);
	  g.setVal(val+distArray.dists[i]->getGradInput(param), 
		   distArray.distIndex[i]);
	}  
    }
  virtual double priorLogProb() const
    {
      double L = 0.0;
      double param=0.0;
      for(int i=0; i<distArray.distIndex.size(); i++)
	{
	  param = getParam(distArray.distIndex[i]);
	  L+=distArray.dists[i]->logProb(param);
	}
      return L;
    }

  // These are non-modifiable methods.
  inline const int getNumPriors() const
    {
      return distArray.getNumDists();
    }
  inline CDist* getPrior(const int ind) const
    {
      assert(ind>=0);
      assert(ind<getNumPriors());
      return distArray.dists[ind];
    }
  inline string getPriorType(const int ind) const
    {
      return distArray.getDistType(ind);
    }
  inline int getPriorIndex(const int ind) const
    {
      return distArray.getDistIndex(ind);
    }
  inline double getPriorGradInput(const double val, const int ind) const
    {
      return distArray.dists[ind]->getGradInput(val);
    }
  void addPrior(CDist* dist, const int index)
    {
      assert(index>=0);
      assert(index<getNumParams());
      distArray.distIndex.push_back(index);
      distArray.dists.push_back(dist);
    }

  void clearPriors()
    {
      distArray.distIndex.clear();
      distArray.dists.clear();
    }

  mxArray* distsToMxArray() const
    {
      return distArray.toMxArray();
    }
  void distsFromMxArray(const mxArray* matlabArray) 
    {
      distArray.fromMxArray(matlabArray);
    }
 private:
  CParamPriors distArray;

};



#endif

    
