#ifndef CDIST_H
#define CDIST_H
#include <cmath>
#include <vector>
#include "CMatrix.h"
#include "CNdlInterfaces.h"
#include "CTransform.h"
#include "ndlutil.h"
              

const string DISTVERSION="0.1";

// Base distribution class.
class CDist : public CMatInterface, public CStreamInterface, public CTransformable {
  
 public:
  CDist(){}
  virtual ~CDist(){}
  unsigned int getNumParams() const
  {
    return nParams;
  }
  void setNumParams(unsigned int num)
  {
    nParams = num;
  }
  virtual double getParam(unsigned int paramNo) const=0;
  virtual void setParam(double val, unsigned int paramNo)=0;
  virtual void getGradParams(CMatrix& g) const
  {
    // This is a dummy function
    cerr << "getGradParams should not be used in CDist" << endl;
    exit(1);
  }
  string getBaseType() const
  {
    return "dist";
  }
  
  
#ifdef _NDLMATLAB
  // returns an mxArray of the dist for use with matlab.
  virtual mxArray* toMxArray() const;
  virtual void fromMxArray(const mxArray* matlabArray);
  // Adds parameters to the mxArray.
  virtual void addParamToMxArray(mxArray* matlabArray) const;
  // Gets the parameters from the mxArray.
  virtual void extractParamFromMxArray(const mxArray* matlabArray);
#endif /* _NDLMATLAB*/
  bool equals(const CDist& dist, double tol=ndlutil::MATCHTOL) const;
 
  virtual void writeParamsToStream(ostream& out) const;
  virtual void readParamsFromStream(istream& in);
  //CDist(CDist& dist);
  // get the gradient with respect to an input.
  virtual double getGradInput(double x) const=0;
  // get the gradient with respect to a matrix of inputs
  virtual void getGradInputs(CMatrix& g, const CMatrix& x)
    {
      DIMENSIONMATCH(g.getRows()==x.getRows());
      DIMENSIONMATCH(g.getCols()==x.getCols());
      for(unsigned int i=0; i<g.getRows(); i++)
	for(unsigned int j=0; j<g.getCols(); j++)
	  g.setVal(getGradInput(x.getVal(i, j)), i, j);
    }
  virtual void setInitParam()=0;
  // Get log probability at a particualar value
  virtual double logProb(double val) const=0;
  virtual double logProb(const CMatrix& x) const
    {
      double ll = 0.0;
      for(unsigned int i=0; i<x.getRows(); i++)
	for(unsigned int j=0; j<x.getCols(); j++)
		ll+=logProb(x.getVal(i, j));
      return ll;
    }
  void setParamName(const string name, unsigned int index)
    {
      
      BOUNDCHECK(index<nParams);
      if(paramNames.size() == index)
	paramNames.push_back(name);
      else 
	{
	  if(paramNames.size()<index)
	    paramNames.resize(index+1, "no name");
	  paramNames[index] = name;
	}
    }
  virtual string getParamName(unsigned int index) const
    {
      
      BOUNDCHECK(index<paramNames.size());
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
  void _init();
  unsigned int nParams;
  string type;
  string distName;
  vector<string> paramNames;
};

// Read and write dists to streams
void writeDistToStream(const CDist& dist, ostream& out);
CDist* readDistFromStream(istream& in);

// The Gaussian distribution.
class CGaussianDist : public CDist {

 public:
  CGaussianDist();
  CGaussianDist(const CGaussianDist&);
  ~CGaussianDist();
  CGaussianDist* clone() const
    {
      return new CGaussianDist(*this);
    }
  double getParam(unsigned int paramNo) const;
  void setParam(double val, unsigned int paramNo);
  double getGradInput(double x) const;
  void setInitParam();
  double logProb(double val) const;

 private:
   void _init();
  double precision;
};

// The gamma distribution
class CGammaDist : public CDist {

 public:
  CGammaDist();
  CGammaDist(const CGammaDist&);
  ~CGammaDist();
  CGammaDist* clone() const
    {
      return new CGammaDist(*this);
    }
  double getParam(unsigned int paramNo) const;
  void setParam(double val, unsigned int paramNo);
  double getGradInput(double x) const;
  void setInitParam();
  double logProb(double val) const;

 private:
   void _init();
  double a;
  double b;
};
// A class which stores distributions in a container for priors over parameters.
// An unusual prior used by Wang in the GPDM thesis.
class CWangDist : public CDist {

 public:
  CWangDist();
  CWangDist(const CWangDist&);
  ~CWangDist();
  CWangDist* clone() const
    {
      return new CWangDist(*this);
    }
  double getParam(unsigned int paramNo) const;
  void setParam(double val, unsigned int paramNo);
  double getGradInput(double x) const;
  void setInitParam();
  double logProb(double val) const;

 private:
   void _init();
  double M;
};
// A class which stores distributions in a container for priors over parameters.
class CParamPriors : CMatInterface {
  
 public:
#ifdef _NDLMATLAB
  mxArray* toMxArray() const;
  void fromMxArray(const mxArray* distArray);
#endif
  void addDist(CDist* dist, unsigned int index)
    {
      
      distIndex.push_back(index);
      dists.push_back(dist);
    }

  void clearDists()
    {
      distIndex.clear();
      dists.clear();
    }
  inline string getDistType(unsigned int ind) const
    {
      
      BOUNDCHECK(ind<getNumDists());
      return dists[ind]->getType();
    }
  inline unsigned int getDistIndex(unsigned int ind) const
    {
      
      BOUNDCHECK(ind<getNumDists());
      return distIndex[ind];
    }
  inline unsigned int getNumDists() const
    {
      return dists.size();
    }
  vector<CDist*> dists;
  vector<int> distIndex;

};

// A virtual base class which makes its descendents regluarisable.
class CRegularisable {

 public:
   virtual ~CRegularisable() {}

  // these are the pure virtual functions.
  virtual unsigned int getNumParams() const=0;
  virtual double getParam(unsigned int paramNo) const=0;
  virtual void setParam(double val, unsigned int paramNo)=0;
  virtual void getGradParams(CMatrix& g) const=0;

  // these are default implementations.
  virtual void getParams(CMatrix& params) const
    {
      DIMENSIONMATCH(params.getRows()==1);
      DIMENSIONMATCH(params.getCols()==getNumParams());
      for(unsigned int i=0; i<params.getCols(); i++)
	params.setVal(getParam(i), i);
    }
  virtual void setParams(const CMatrix& params)
    {
      DIMENSIONMATCH(params.getRows()==1);
      DIMENSIONMATCH(params.getCols()==getNumParams());
      for(unsigned int i=0; i<params.getCols(); i++)
	setParam(params.getVal(i), i);
    }
  
  virtual void addPriorGrad(CMatrix& g) const
    {
      DIMENSIONMATCH(g.getRows()==1);
      DIMENSIONMATCH(g.getCols()==getNumParams());
      double param=0.0;
      for(unsigned int i=0; i<distArray.distIndex.size(); i++)
	{
	  param=getParam(distArray.distIndex[i]);
	  g.addVal(getPriorGradInput(param, i), 
		   distArray.distIndex[i]);
	}  
    }
      
  virtual void writePriorsToStream(ostream& out) const
    {
      for(unsigned int i=0; i<distArray.distIndex.size(); i++)
	{
	  out << "priorIndex=" << distArray.distIndex[i] << endl;
	  writeDistToStream(*distArray.dists[i], out);
	}
    }
  virtual void readPriorsFromStream(istream& in, unsigned int numPriors)
    {
      string line;
      vector<string> tokens;
      for(unsigned int i=0; i<numPriors; i++)
	{
	  CDist* prior;
	  getline(in, line);
	  ndlstrutil::tokenise(tokens, line, "=");
	  if(tokens.size()>2 || tokens[0]!="priorIndex")
	    throw ndlexceptions::StreamFormatError("priorIndex");
	  prior = readDistFromStream(in);
	  addPrior(prior, atol(tokens[1].c_str()));
	  tokens.clear();
	}
    }

  virtual double priorLogProb() const
    {
      double L = 0.0;
      double param=0.0;
      for(unsigned int i=0; i<distArray.distIndex.size(); i++)
	{
	  param = getParam(distArray.distIndex[i]);
	  L+=distArray.dists[i]->logProb(param);
	}
      return L;
    }

  // These are non-modifiable methods.
  inline unsigned int getNumPriors() const
    {
      return distArray.getNumDists();
    }
  inline CDist* getPrior(unsigned int ind) const
    {
      
      BOUNDCHECK(ind<getNumPriors());
      return distArray.dists[ind];
    }
  inline string getPriorType(unsigned int ind) const
    {
      return distArray.getDistType(ind);
    }
  inline unsigned int getPriorIndex(unsigned int ind) const
    {
      return distArray.getDistIndex(ind);
    }
  inline double getPriorGradInput(double val, unsigned int ind) const
    {
      return distArray.dists[ind]->getGradInput(val);
    }
  void addPrior(CDist* dist, unsigned int index)
    {
      
      BOUNDCHECK(index<getNumParams());
      distArray.distIndex.push_back(index);
      distArray.dists.push_back(dist);
    }

  void clearPriors()
    {
      distArray.distIndex.clear();
      distArray.dists.clear();
    }
#ifdef _NDLMATLAB
  mxArray* distsToMxArray() const
    {
      return distArray.toMxArray();
    }
  void distsFromMxArray(const mxArray* matlabArray) 
    {
      distArray.fromMxArray(matlabArray);
    }
#endif
 private:
  CParamPriors distArray;

};



#endif

    
