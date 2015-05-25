#ifndef CKERN_H
#define CKERN_H
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "ndlassert.h"
#include "CTransform.h"
#include "CDataModel.h"
#include "CMatrix.h"
#include "CDist.h"
#include "ndlstrutil.h"
#include "ndlutil.h"
using namespace std;


const string KERNVERSION="0.1";

// The base class for all kernels. It implements: CMatinterface, for
// saving and loading from MATLAB; CRegularisable, for placing priors
// over parameters; CTransformable, for allowing parameters to be
// transformed so that they are only optimised in, for example,
// positive half-spaces

class CKern : public CMatInterface, public CStreamInterface, public CTransformable, public CRegularisable {
 public:
  CKern() : updateXused(false) {}
  CKern(const CMatrix& X) : updateXused(false) {}
  CKern(unsigned int inDim) : updateXused(false) {}
  CKern(const CKern& kern) : updateXused(false) {}
  
  virtual ~CKern()  {}
  virtual CKern* clone() const=0;
  // set initial parameters.
  virtual void setInitParam()=0;
  // compute an element of the diagonal.
  virtual double diagComputeElement(const CMatrix& X, unsigned int index) const=0;
  // compute the entire diagonal
  virtual void diagCompute(CMatrix& d, const CMatrix& X) const
  {
    DIMENSIONMATCH(X.rowsMatch(d));
    DIMENSIONMATCH(d.getCols()==1);
    for(unsigned int i=0; i<X.getRows(); i++)
      d.setVal(diagComputeElement(X, i), i);
  }
  // Compute the diagonal at particular indices.
  virtual void diagCompute(CMatrix& d, const CMatrix& X, const vector<unsigned int> indices) const
  {
    DIMENSIONMATCH(d.getRows()==indices.size());
    DIMENSIONMATCH(d.getCols()==1);
    for(unsigned int i=0; i<indices.size(); i++)
      d.setVal(diagComputeElement(X, indices[i]), i);
  }
 
  // Set the parameters of the kernel.
  virtual void setParam(double, unsigned int)=0;
  // Get gradients of the kernel with respect to input values.
  //   g[i].val(k,j) = d kern(X_row_i,X2_row_k)/ d x_component_j
  virtual void getGradX(vector<CMatrix*>& gX, const CMatrix& X, const CMatrix& X2, bool addG=false) const
  {
    for(unsigned int i=0; i<X.getRows(); i++)
      getGradX(*gX[i], X, i, X2, addG);
  }
  virtual void getGradX(CMatrix& g, const CMatrix& X, unsigned int pointNo, const CMatrix& X2, bool addG=false) const=0;
  // Get gradients of the kernel diagonal with respect to input values.
  // WVB: I think this is separate from getGradX just because the diagonal vals
  // (i.e. kern(X_row_i,X_row_i) should always be taken from the same matrix X.
  // Not to imply that I understand why we need X and X2 in getGradX, really.
  virtual void getDiagGradX(CMatrix& g, const CMatrix& X, bool addG=false) const=0;
  // return the `signal strength' of the kernel.
  virtual double getVariance() const=0;
  // set the `signal strength' of the kernel.
  virtual void setVariance(double val)=0;
  // Return the white noise component of the kernel.
  virtual double getWhite() const
  {
    return 0.0;
  }
  // Compute an element of the kernel matrix.
  virtual double computeElement(const CMatrix& X1, unsigned int index1,
			 const CMatrix& X2, unsigned int index2) const=0;
  // Compute specified rows and columns of the kernel matrix.
  virtual void compute(CMatrix& K, const CMatrix& X1, const vector<unsigned int> indices1,
			  const CMatrix& X2, const vector<unsigned int> indices2) const
  {
    DIMENSIONMATCH(K.getRows()==indices1.size());
    DIMENSIONMATCH(K.getCols()==indices2.size());
    for(unsigned int i=0; i<indices1.size(); i++)
    {
      for(unsigned int j=0; j<indices2.size(); j++)
      {
	BOUNDCHECK(indices1[i]<X1.getRows());
	BOUNDCHECK(indices2[j]<X2.getRows());
	K.setVal(computeElement(X1, indices1[i], X2, indices2[j]), i, j);
      }
    }
  }
  virtual void compute(CMatrix& K, const CMatrix& X, const vector<unsigned int> indices) const
  {
    DIMENSIONMATCH(K.getRows()==indices.size());
    MATRIXPROPERTIES(K.isSquare());
    double k = 0.0;
    for(unsigned int i=0; i<indices.size(); i++)
    {
      for(unsigned int j=0; j<i; j++)
      {
	BOUNDCHECK(indices[i]<X.getRows());
	BOUNDCHECK(indices[j]<X.getRows());
	k = computeElement(X, indices[i], X, indices[j]);
	K.setVal(k, i, j);
	K.setVal(k, j, i);
      }
      K.setVal(diagComputeElement(X, indices[i]), i, i);
    }
  }
  // Compute the kernel matrix for a data set.
  virtual void compute(CMatrix& K, const CMatrix& X) const
  {
    DIMENSIONMATCH(K.rowsMatch(X));
    MATRIXPROPERTIES(K.isSquare());
    double k = 0.0;
    for(unsigned int i=0; i<K.getRows(); i++)
    {
      for(unsigned int j=0; j<i; j++)
      {
	k = computeElement(X, i, X, j);
	K.setVal(k, i, j);
	K.setVal(k, j, i);
      }
      K.setVal(diagComputeElement(X, i), i, i);
    }	      
    K.setSymmetric(true);
  }
  // Compute portions of the kernel matrix.
  virtual void compute(CMatrix& K, const CMatrix& X, const CMatrix& X2) const
  {
    DIMENSIONMATCH(K.rowsMatch(X));
    DIMENSIONMATCH(K.getCols()==X2.getRows());
    for(unsigned int i=0; i<K.getRows(); i++)
    {
      for(unsigned int j=0; j<K.getCols(); j++)
      {
	K.setVal(computeElement(X, i, X2, j), i, j);
      }
    }	      
  }
  // Compute portions of the kernel matrix.
  virtual void compute(CMatrix& K, const CMatrix& X, const CMatrix& X2, unsigned int row) const
  {
    DIMENSIONMATCH(K.rowsMatch(X));
    DIMENSIONMATCH(K.getCols()==1);
    for(unsigned int i=0; i<K.getRows(); i++)
    {
      K.setVal(computeElement(X, i, X2, row), i, 0);
    }
  }
  // Dummy function to allow CTransformable to be used.
  virtual void getGradParams(CMatrix& g) const
  {
    // This is a dummy function
    cerr << "getGradParams should not be used in CKern" << endl;
    exit(1);
  }
  // Compute the gradient of the kernel matrix with respect to parameters given an additional gradient matrix.
  virtual void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd, bool regularise=true) const
  {
    DIMENSIONMATCH(g.getRows()==1);
    DIMENSIONMATCH(g.getCols()==nParams);
    DIMENSIONMATCH(X.getRows()==cvGrd.getRows());
    DIMENSIONMATCH(X2.getRows()==cvGrd.getCols());
    for(unsigned int i=0; i<nParams; i++)
      g.setVal(getGradParam(i, X, X2, cvGrd), i);
    if(regularise)
      addPriorGrad(g);
  }
  virtual void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd, bool regularise=true) const
  {
    DIMENSIONMATCH(g.getRows()==1);
    DIMENSIONMATCH(g.getCols()==nParams);
    DIMENSIONMATCH(X.rowsMatch(cvGrd));
    DIMENSIONMATCH(cvGrd.isSquare());
    for(unsigned int i=0; i<nParams; i++)
      g.setVal(getGradParam(i, X, cvGrd), i);
    if(regularise)
      addPriorGrad(g); /// don't forget to add prior gradient at the end.
  }
  virtual void getDiagGradParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrad, bool regularise=true) const
  {
    //TODO Code explicitly for things like RBF where it will always be zero.
    CMatrix xi(1, X.getCols());
    CMatrix cvGradi(1, 1);
    CMatrix gtemp(1, g.getCols());
    g.zeros();
    for(unsigned int i=0; i<X.getRows(); i++)
    {
      xi.copyRowRow(0, X, i);
      cvGradi.copyRowRow(0, cvGrad, i);
      cvGradi.setSymmetric(true);
      getGradParams(gtemp, xi, cvGradi, regularise);
      g.axpy(gtemp, 1.0);
    }
  }
  // Get gradient of a particular parameter.
  virtual double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd) const=0;
  virtual double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& cvGrd) const=0;
  // Called to indicate the value of X has changed and kernel should do any
  // precomputation it needs to do per value of X.
  virtual void updateX(const CMatrix& X) {}

  virtual bool isUpdateXused() const
  {
    return updateXused;
  }
  virtual bool setUpdateXused(bool val) 
  {
    updateXused = val;
  }
  // For compound kernels when a new kernel is added.
  virtual unsigned int addKern(const CKern* kern)
  {
    cerr << "You cannot add a kernel to this kernel." << endl;
    return 0;
  }
  // Get a particular parameter.
  virtual double getParam(unsigned int) const=0;
  // Test if kernel leads to stationary functions.
  bool isStationary() const
  {
    return stationary;
  }
  void setStationary(bool val)
  {
    stationary = val;
  }
  // Set the parameters from a vector of parameters.
  void setParams(const CMatrix& paramVec)
    {
      for(unsigned int i=0; i<nParams; i++)
	setParam(paramVec.getVal(i), i);
    }
  // Place the parameters in a vector.
  void getParams(CMatrix& paramVec) const
    {
      for(unsigned int i=0; i<nParams; i++)
	paramVec.setVal(getParam(i), i);
    }
  // Return a string representing the kernel type.
  inline string getType() const
  {
      return type;
  }
  // Set a string representing the kernel type.
  inline void setType(const string name)
  {
    type = name;
  }
  string getBaseType() const
  {
    return "kern";
  }
  // Get the long name of the kernel.
  inline string getName() const
    {
      return kernName;
    }
  // Set the long name of the kernel.
  inline void setName(const string name)
    {
      kernName = name;
    }
  // Set the input dimension.
  inline void setInputDim(unsigned int dim)
    {
      inputDim = dim;
      setInitParam();
    }
  // Get the input dimension.
  inline unsigned getInputDim() const
    {
      return inputDim;
    }
  // How many kernel parameters are there?
  inline unsigned int getNumParams() const
  {
    return nParams;
  }
  // Assign a name to the kernel parameters.
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
  // Get the name of a kernel parameter.
  virtual string getParamName(unsigned int index) const
    {
      BOUNDCHECK(index<paramNames.size());
      return paramNames[index];
    }
  // Write out the kernel parameters to a stream.
  virtual void writeParamsToStream(ostream& out) const;
  // Read in the kernel parameters from a stream.
  virtual void readParamsFromStream(istream& in);
  // Get the gradients of the parameters associated with the priors.
  void getGradPrior(CMatrix& g) const;
  // Get the log probabilities associated with the priors.
  void getPriorLogProb(CMatrix& L) const;
  // Display the kernel on an ostream.
  virtual ostream& display(ostream& os) const;
  
#ifdef _NDLMATLAB
  // Create a kernel from an mxArray* object
  CKern(mxArray* kern){}
  // returns an mxArray of the kern for use with matlab.
  virtual mxArray* toMxArray() const;
  virtual void fromMxArray(const mxArray* matlabArray);
  // Adds parameters to the mxArray.
  virtual void addParamToMxArray(mxArray* matlabArray) const;
  // Gets the parameters from the mxArray.
  virtual void extractParamFromMxArray(const mxArray* matlabArray);
#endif /* _NDLMATLAB*/
  // Get the gradient of the transformed parameters.
  void getGradTransParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd, bool regularise=true) const;
  void getGradTransParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd, bool regularise=true) const;
  void getDiagGradTransParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd, bool regularise=true) const;
  // specify tests for equality between kernels.
  bool equals(const CKern& kern, double tol=ndlutil::MATCHTOL) const;

 protected:
  unsigned int nParams;
  string kernName;
  string type;
  vector<string> paramNames;
 private:
  bool updateXused;
  unsigned int inputDim;
  bool stationary;
};

// CArdKern is the base class for any kernel that uses multiple input parameters.
class CArdKern : public CKern {
 public:
  CArdKern() : CKern() {}
  CArdKern(const CMatrix& X) : CKern(X) {}
  CArdKern(unsigned int inDim) : CKern(inDim) {}
  CArdKern(const CKern& kern) : CKern(kern) {}
#ifdef _NDLMATLAB
  virtual void addParamToMxArray(mxArray* matlabArray) const;
  // Gets the parameters from the mxArray.
  virtual void extractParamFromMxArray(const mxArray* matlabArray);
  // returns sum(sum(cvGrd.*dK/dparam)) 
#endif
 protected: 
  CMatrix scales;
};

// Component kernel (such as cmpnd or tensor)
class CComponentKern : public CKern 
{
 public:
  CComponentKern() : CKern() {}
  CComponentKern(unsigned int inDim) : CKern(inDim) {}
  CComponentKern(const CMatrix& X) : CKern(X) {}
  CComponentKern(const CComponentKern& kern) : CKern(kern), components(kern.components) {}
  virtual unsigned int addKern(const CKern* kern)
  {
    components.push_back(kern->clone());
    unsigned int oldNParams = nParams;
    nParams+=kern->getNumParams();
    for(size_t i=0; i<kern->getNumTransforms(); i++)
      addTransform(kern->getTransform(i), kern->getTransformIndex(i)+oldNParams);      
    setStationary(isStationary() && kern->isStationary());
    return components.size()-1;
  }
  virtual void setParam(double val, unsigned int paramNo)
  {
    unsigned int start = 0;
    unsigned int end = 0;
    for(size_t i=0; i<components.size(); i++)
    {
      end = start+components[i]->getNumParams()-1;
      if(paramNo <= end)
      {
	components[i]->setParam(val, paramNo-start);
	return;
      }      
      start = end + 1;
    }
  }
  // Parameters are kernel parameters
  virtual double getParam(unsigned int paramNo) const
  {
    unsigned int start = 0;
    unsigned int end = 0;
    for(size_t i=0; i<components.size(); i++)
    {
      end = start+components[i]->getNumParams()-1;
      if(paramNo <= end)
	return components[i]->getParam(paramNo-start);
      start = end + 1;
    }
    return -1;
  }
  virtual string getParamName(unsigned int paramNo) const
  {
    unsigned int start = 0;
    unsigned int end = 0;
    for(size_t i=0; i<components.size(); i++)
    {
      end = start+components[i]->getNumParams()-1;
      if(paramNo <= end)
	return components[i]->getType() + components[i]->getParamName(paramNo-start);
      start = end + 1;
    }
    return "";
  }
  virtual void updateX(const CMatrix& X)
  {
    for(size_t i=0; i<components.size(); i++)
      components[i]->updateX(X);
  }
  virtual void addPrior(CDist* prior, unsigned int index) 
  {
    throw ndlexceptions::Error("Error cannot add priors to component kernels directly, please add to the components.");
  }

  virtual double priorLogProb() const
  {
    double L = 0.0;
    for(unsigned int i=0; i<components.size(); i++)
    {
      L+=components[i]->priorLogProb();
    }
    return L;
  }
  
  virtual void readParamsFromStream(istream& in); 
  virtual void writeParamsToStream(ostream& out) const;
  virtual unsigned int getNumKerns() const
  {
    return components.size();
  }

  
  
#ifdef _NDLMATLAB
  // sets the parameters in the mxArray.
  virtual void addParamToMxArray(mxArray* matlabArray) const;
  // Gets the parameters from the mxArray.
  virtual void extractParamFromMxArray(const mxArray* matlabArray); 
#endif
 protected:
  // this is a heterogeneous container.
  vector<CKern*> components;
  
};
// Compound Kernel --- This kernel combines other kernels additively together.
class CCmpndKern: public CComponentKern {
 public:
  CCmpndKern();
  CCmpndKern(unsigned int inDim);
  CCmpndKern(const CMatrix& X);
  ~CCmpndKern();
  CCmpndKern(const CCmpndKern&);
  CCmpndKern* clone() const
  {
    return new CCmpndKern(*this);
  }
  //CCmpndKern(vector<CKern*> kernels);
  
  double getVariance() const;
  void setVariance(double val)
  {
    double totalVariance = getVariance();
    double factor = val/totalVariance;
    for(size_t i=0; i<components.size(); i++)
    {
      double newVariance = components[i]->getVariance()*factor;
      components[i]->setVariance(newVariance);
    }  
  }
  double getWhite() const;
  void setInitParam();
  double diagComputeElement(const CMatrix& X, unsigned int index1) const;
  void diagCompute(CMatrix& d, const CMatrix& X) const;
  void getGradX(CMatrix& g, const CMatrix& X, unsigned int pointNo, const CMatrix& X2, bool addG=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, bool addG=false) const;
  double computeElement(const CMatrix& X1, unsigned int index1, 
		 const CMatrix& X2, unsigned int index2) const;
/*   void compute(CMatrix& K, const CMatrix& X) const; */
/*   void compute(CMatrix& K, const CMatrix& X, const CMatrix& X2) const; */
/*   void compute(CMatrix& K, const CMatrix& X, const CMatrix& X2, unsigned int row) const; */
  //  void compute(CMatrix& K, const CMatrix& X, const vector<unsigned int> indices) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd, bool regularise=true) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd, bool regularise=true) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& cvGrd) const;
private:
  void _init();
};

// Tensor Kernel --- This kernel combines other multiplicitavely together.
class CTensorKern: public CComponentKern {
 public:
  CTensorKern();
  CTensorKern(unsigned int inDim);
  CTensorKern(const CMatrix& X);
  ~CTensorKern();
  CTensorKern(const CTensorKern&);
  CTensorKern* clone() const
  {
    return new CTensorKern(*this);
  }
  CTensorKern(const CTensorKern&, unsigned int i);
  //CTensorKern(vector<CKern*> kernels);
  
  double getVariance() const;
  void setVariance(double val)
  {
    double totalVariance = getVariance();
    double factor = val/totalVariance;
    for(size_t i=0; i<components.size(); i++)
    {
      double newVariance = components[i]->getVariance()*factor;
      components[i]->setVariance(newVariance);
    } 
  }
  double getWhite() const;
  void setInitParam();
  double diagComputeElement(const CMatrix& X, unsigned int index1) const;
  void diagCompute(CMatrix& d, const CMatrix& X) const;
  void getGradX(CMatrix& g, const CMatrix& X, unsigned int pointNo, const CMatrix& X2, bool addG=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, bool addG=false) const;
  double computeElement(const CMatrix& X1, unsigned int index1, const CMatrix& X2, unsigned int index2) const;
/*   void compute(CMatrix& K, const CMatrix& X) const; */
/*   void compute(CMatrix& K, const CMatrix& X, const CMatrix& X2) const; */
/*   void compute(CMatrix& K, const CMatrix& X, const CMatrix& X2, unsigned int row) const; */
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd, bool regularise=true) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd, bool regularise=true) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& cvGrd) const;
  unsigned int addKern(const CKern* kern);
private:
  void _init();
};

// White Noise Kernel.
class CWhiteKern: public CKern {
 public:
  CWhiteKern();
  CWhiteKern(unsigned int inDim);
  CWhiteKern(const CMatrix& X);
  ~CWhiteKern();
  CWhiteKern(const CWhiteKern&);
  CWhiteKern* clone() const
  {
    return new CWhiteKern(*this);
  }
  double getVariance() const;
  void setVariance(double val)
  {
    variance = val;
  }
  void setInitParam();
  double diagComputeElement(const CMatrix& X, unsigned int index) const;
  void diagCompute(CMatrix& d, const CMatrix& X) const;
  void setParam(double val, unsigned int paramNum);
  double getParam(unsigned int paramNum) const;
  void getGradX(CMatrix& g, const CMatrix& X, unsigned int pointNo, const CMatrix& X2, bool addG=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, bool addG=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, unsigned int index1, 
			const CMatrix& X2, unsigned int index2) const;
  void compute(CMatrix& K, const CMatrix& X) const;
  void compute(CMatrix& K, const CMatrix& X, const CMatrix& X2) const;
  void compute(CMatrix& K, const CMatrix& X, const CMatrix& X2, unsigned int row) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& cvGrd) const;


 private:
  void _init();
  double variance;
};
// Whitefixed Noise Kernel.
class CWhitefixedKern: public CKern {
 public:
  CWhitefixedKern();
  CWhitefixedKern(unsigned int inDim);
  CWhitefixedKern(const CMatrix& X);
  ~CWhitefixedKern();
  CWhitefixedKern(const CWhitefixedKern&);
  CWhitefixedKern* clone() const
  {
    return new CWhitefixedKern(*this);
  }
  double getVariance() const;
  void setVariance(double val)
  {
    variance = val;
  }
  void setInitParam();
  double diagComputeElement(const CMatrix& X, unsigned int index) const;
  void diagCompute(CMatrix& d, const CMatrix& X) const;
  void setParam(double val, unsigned int paramNum);
  double getParam(unsigned int paramNum) const;
  void getGradX(CMatrix& g, const CMatrix& X, unsigned int pointNo, const CMatrix& X2, bool addG=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, bool addG=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, unsigned int index1, 
			const CMatrix& X2, unsigned int index2) const;
  void compute(CMatrix& K, const CMatrix& X) const;
  void compute(CMatrix& K, const CMatrix& X, const CMatrix& X2) const;
  void compute(CMatrix& K, const CMatrix& X, const CMatrix& X2, unsigned int row) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& cvGrd) const;
  void writeParamsToStream(ostream& out) const;
  void readParamsFromStream(istream& in);
  ostream& display(ostream& os) const;


 private:
  void _init();
  double variance;
};

// Bias Kernel.
class CBiasKern: public CKern {
 public:
  CBiasKern();
  CBiasKern(unsigned int inDim);
  CBiasKern(const CMatrix& X);
  ~CBiasKern();
  CBiasKern(const CBiasKern&);
  CBiasKern* clone() const
  {
    return new CBiasKern(*this);
  }
  double getVariance() const;
  void setVariance(double val)
  {
    variance = val;
  }
  void setInitParam();
  double diagComputeElement(const CMatrix& X, unsigned int index) const;
  void diagCompute(CMatrix& d, const CMatrix& X) const;
  void setParam(double val, unsigned int paramNum);
  double getParam(unsigned int paramNum) const;
  void getGradX(CMatrix& g, const CMatrix& X, unsigned int pointNo, const CMatrix& X2, bool addG=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, bool addG=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, unsigned int index1,
		 const CMatrix& X2, unsigned int index2) const;
  void compute(CMatrix& K, const CMatrix& X) const;
  void compute(CMatrix& K, const CMatrix& X, const CMatrix& X2) const;
  void compute(CMatrix& K, const CMatrix& X, const CMatrix& X2, unsigned int row) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& cvGrd) const;



 private:
  void _init();
  double variance;

};
// RBF Kernel, also known as the Gaussian or squared exponential kernel.
class CRbfKern: public CKern {
 public:
  CRbfKern();
  CRbfKern(unsigned int inDim);
  CRbfKern(const CMatrix& X);
  ~CRbfKern();
  CRbfKern(const CRbfKern&);
  CRbfKern* clone() const
  {
    return new CRbfKern(*this);
  }
  double getVariance() const;
  void setVariance(double val)
  {
    variance = val;
  }
  void setInverseWidth(double val)
  {
    inverseWidth = val;
  }
  double getInverseWidth() const
  {
    return inverseWidth;
  }
  void setLengthScale(double val)
  {
    inverseWidth = 1/(val*val);
  }
  double getLengthScale() const
  {
    return 1/sqrt(inverseWidth);
  }
  void setInitParam();
  double diagComputeElement(const CMatrix& X, unsigned int index) const;
  void diagCompute(CMatrix& d, const CMatrix& X) const;
  void setParam(double val, unsigned int paramNum);
  double getParam(unsigned int paramNum) const;
  void getGradX(CMatrix& g, const CMatrix& X, unsigned int pointNo, const CMatrix& X2, bool addG=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, bool addG=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, unsigned int index1, 
			const CMatrix& X2, unsigned int index2) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd, bool regularise=true) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd, bool regularise=true) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& cvGrd) const;
  void updateX(const CMatrix& X);
  
 private:
  void _init();
  double variance;
  double inverseWidth;
  mutable CMatrix Xdists;
  
};

// Rational Quadratic Kernel
class CRatQuadKern: public CKern {
 public:
  CRatQuadKern();
  CRatQuadKern(unsigned int inDim);
  CRatQuadKern(const CMatrix& X);
  ~CRatQuadKern();
  CRatQuadKern(const CRatQuadKern&);
  CRatQuadKern* clone() const
  {
    return new CRatQuadKern(*this);
  }
  double getVariance() const;
  void setVariance(double val)
  {
    variance = val;
  }
  void setInitParam();
  double diagComputeElement(const CMatrix& X, unsigned int index) const;
  void diagCompute(CMatrix& d, const CMatrix& X) const;
  void setParam(double val, unsigned int paramNum);
  double getParam(unsigned int paramNum) const;
  void getGradX(CMatrix& g, const CMatrix& X, unsigned pointNo, const CMatrix& X2, bool addG=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, bool addG=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, unsigned int index1, 
			const CMatrix& X2, unsigned int index2) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd, bool regularise=true) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd, bool regularise=true) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& cvGrd) const;
  void updateX(const CMatrix& X);
  
 private:
  void _init();
  double variance;
  double alpha;
  double lengthScale;
  mutable CMatrix Xdists;
  
};

// Matern kernel with dof=3/2.
class CMatern32Kern: public CKern {
 public:
  CMatern32Kern();
  CMatern32Kern(unsigned int inDim);
  CMatern32Kern(const CMatrix& X);
  ~CMatern32Kern();
  CMatern32Kern(const CMatern32Kern&);
  CMatern32Kern* clone() const
  {
    return new CMatern32Kern(*this);
  }
  double getVariance() const;
  void setVariance(double val)
  {
    variance = val;
  }
  void setInitParam();
  double diagComputeElement(const CMatrix& X, unsigned int index) const;
  void diagCompute(CMatrix& d, const CMatrix& X) const;
  void setParam(double val, unsigned int paramNum);
  double getParam(unsigned int paramNum) const;
  void getGradX(CMatrix& g, const CMatrix& X, unsigned int pointNo, const CMatrix& X2, bool addG=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, bool addG=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, unsigned int index1, 
			const CMatrix& X2, unsigned int index2) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd, bool regularise=true) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd, bool regularise=true) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& cvGrd) const;
  void updateX(const CMatrix& X);
  
 private:
  void _init();
  double variance;
  double lengthScale;
  mutable CMatrix Xdists;
  
};

// Matern kernel with dof=5/2.
class CMatern52Kern: public CKern {
 public:
  CMatern52Kern();
  CMatern52Kern(unsigned int inDim);
  CMatern52Kern(const CMatrix& X);
  ~CMatern52Kern();
  CMatern52Kern(const CMatern52Kern&);
  CMatern52Kern* clone() const
  {
    return new CMatern52Kern(*this);
  }
  double getVariance() const;
  void setVariance(double val)
  {
    variance = val;
  }
  void setInitParam();
  double diagComputeElement(const CMatrix& X, unsigned int index) const;
  void diagCompute(CMatrix& d, const CMatrix& X) const;
  void setParam(double val, unsigned int paramNum);
  double getParam(unsigned int paramNum) const;
  void getGradX(CMatrix& g, const CMatrix& X, unsigned int pointNo, const CMatrix& X2, bool addG=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, bool addG=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, unsigned int index1, 
			const CMatrix& X2, unsigned int index2) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd, bool regularise=true) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd, bool regularise=true) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& cvGrd) const;
  void updateX(const CMatrix& X);
  
 private:
  void _init();
  double variance;
  double lengthScale;
  mutable CMatrix Xdists;
  
};

// Linear Kernel, also known as the inner product kernel.
class CLinKern: public CKern {
 public:
  CLinKern();
  CLinKern(unsigned int inDim);
  CLinKern(const CMatrix& X);
  ~CLinKern();
  CLinKern(const CLinKern&);
  CLinKern* clone() const
  {
    return new CLinKern(*this);
  }
  double getVariance() const;
  void setVariance(double val)
  {
    variance = val;
  }
  void setInitParam();
  double diagComputeElement(const CMatrix& X, unsigned int index) const;
  void setParam(double val, unsigned int paramNum);
  double getParam(unsigned int paramNum) const;
  void getGradX(CMatrix& g, const CMatrix& X, unsigned int pointNo, const CMatrix& X2, bool addG=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, bool addG=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, unsigned int index1, 
			const CMatrix& X2, unsigned int index2) const;
  void compute(CMatrix& K, const CMatrix& X) const;
  void compute(CMatrix& K, const CMatrix& X, const CMatrix& X2) const;
  void compute(CMatrix& K, const CMatrix& X, const CMatrix& X2, unsigned int row) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& cvGrd) const;


 private:
  void _init();
  double variance;
};

// MLP Kernel or arcsin kernel. Based on a multi-layer perceptron with infinite hidden nodes. See Williams (1996) "Computing with Infinite Networks" in NIPS 9.
class CMlpKern: public CKern {
 public:
  CMlpKern();
  CMlpKern(unsigned int inDim);
  CMlpKern(const CMatrix& X);
  ~CMlpKern();
  CMlpKern(const CMlpKern&);
  CMlpKern* clone() const
  {
    return new CMlpKern(*this);
  }
  double getVariance() const;
  void setVariance(double val)
  {
    variance = val;
  }
  void setInitParam();
  double diagComputeElement(const CMatrix& X, unsigned int index) const;
  void setParam(double val, unsigned int paramNum);
  double getParam(unsigned int paramNum) const;
  void getGradX(CMatrix& g, const CMatrix& X, unsigned int pointNo, const CMatrix& X2, bool addG=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, bool addG=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, unsigned int index1, 
			const CMatrix& X2, unsigned int index2) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd, bool regularise=true) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd, bool regularise=true) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& cvGrd) const;
  
 private:
  void _init();
  double weightVariance;
  double biasVariance;
  double variance;
  mutable CMatrix innerProdi;
  mutable CMatrix innerProdj;
};

// Polynomial Kernel, not generally recommended as it `extreme behaviour' outside the region where the argument's absolute value is less than 1.
class CPolyKern: public CKern {
 public:
  CPolyKern();
  CPolyKern(unsigned int inDim);
  CPolyKern(const CMatrix& X);
  ~CPolyKern();
  CPolyKern(const CPolyKern&);
  CPolyKern* clone() const
  {
    return new CPolyKern(*this);
  }
  void setDegree(double val)
  {
    degree = val;
  }
  double getDegree() const
  {
    return degree;
  }
  double getVariance() const;
  void setVariance(double val)
  {
    variance = val;
  }
  void setInitParam();
  double diagComputeElement(const CMatrix& X, unsigned int index) const;
  void setParam(double val, unsigned int paramNum);
  double getParam(unsigned int paramNum) const;
  void getGradX(CMatrix& g, const CMatrix& X, unsigned int pointNo, const CMatrix& X2, bool addG=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, bool addG=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, unsigned int index1, 
			const CMatrix& X2, unsigned int index2) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd, bool regularise=true) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd, bool regularise=true) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& cvGrd) const;
  void writeParamsToStream(ostream& out) const;
  void readParamsFromStream(istream& in);
#ifdef _NDLMATLAB
  void addParamToMxArray(mxArray* matlabArray) const;
  void extractParamFromMxArray(const mxArray* matlabArray);
#endif
 private:
  void _init();
  double weightVariance;
  double biasVariance;
  double variance;
  double degree;
  mutable CMatrix innerProdi;
};

// Linear ARD Kernel --- automatic relevance determination version of the linear kernel.
class CLinardKern: public CArdKern {
 public:
  CLinardKern();
  CLinardKern(unsigned int inDim);
  CLinardKern(const CMatrix& X);
  ~CLinardKern();
  CLinardKern(const CLinardKern&);
  CLinardKern* clone() const
  {
    return new CLinardKern(*this);
  }
  double getVariance() const;
  void setVariance(double val)
  {
    variance = val;
  }
  void setInitParam();
  double diagComputeElement(const CMatrix& X, unsigned int index) const;
  void setParam(double val, unsigned int paramNum);
  double getParam(unsigned int paramNum) const;
  void getGradX(CMatrix& g, const CMatrix& X, unsigned int pointNo, const CMatrix& X2, bool addG=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, bool addGrad=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, unsigned int index1, 
			const CMatrix& X2, unsigned int index2) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd, bool regularise = true) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd, bool regularise = true) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& cvGrd) const;
  
  
 private:
  void _init();
  double variance;
};

// RBF ARD Kernel --- automatic relevance determination of the RBF kernel.
class CRbfardKern: public CArdKern {
 public:
  CRbfardKern();
  CRbfardKern(unsigned int inDim);
  CRbfardKern(const CMatrix& X);
  ~CRbfardKern();
  CRbfardKern(const CRbfardKern&);
  CRbfardKern* clone() const
  {
    return new CRbfardKern(*this);
  }
  double getVariance() const;
  void setVariance(double val)
  {
    variance = val;
  }
  void setInitParam();
  double diagComputeElement(const CMatrix& X, unsigned int index) const;
  void setParam(double val, unsigned int paramNum);
  double getParam(unsigned int paramNum) const;
  void getGradX(CMatrix& g, const CMatrix& X, unsigned int pointNo, const CMatrix& X2, bool addG=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, bool addGrad=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, unsigned int index1, 
			const CMatrix& X2, unsigned int index2) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd, bool regularise=true) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd, bool regularise=true) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& cvGrd) const;
  
  
 private:
  void _init();
  double variance;
  double inverseWidth;
  mutable CMatrix gscales;
};

// MLP ARD Kernel --- automatic relevance determination version of the MLP kernel.
class CMlpardKern: public CArdKern {
 public:
  CMlpardKern();
  CMlpardKern(unsigned int inDim);
  CMlpardKern(const CMatrix& X);
  ~CMlpardKern();
  CMlpardKern(const CMlpardKern&);
  CMlpardKern* clone() const
  {
    return new CMlpardKern(*this);
  }
  double getVariance() const;
  void setVariance(double val)
  {
    variance = val;
  }
  void setInitParam();
  double diagComputeElement(const CMatrix& X, unsigned int index) const;
  void setParam(double val, unsigned int paramNum);
  double getParam(unsigned int paramNum) const;
  void getGradX(CMatrix& g, const CMatrix& X, unsigned int pointNo, const CMatrix& X2, bool addG=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, bool addGrad=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, unsigned int index1, 
			const CMatrix& X2, unsigned int index2) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd, bool regularise=true) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd, bool regularise=true) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& cvGrd) const;
  
  
 private:
  void _init();
  double weightVariance;
  double biasVariance;
  double variance;
  mutable CMatrix innerProdi;
  mutable CMatrix innerProdj;
  mutable CMatrix gscales;
};

// Polynomial ARD Kernel --- automatic relevance determination version of the polynomial kernel.
class CPolyardKern: public CArdKern {
 public:
  CPolyardKern();
  CPolyardKern(unsigned int inDim);
  CPolyardKern(const CMatrix& X);
  ~CPolyardKern();
  CPolyardKern(const CPolyardKern&);
  CPolyardKern* clone() const
  {
    return new CPolyardKern(*this);
  }
  void setDegree(double val)
  {
    degree = val;
  }
  double getDegree() const
  {
    return degree;
  }
  double getVariance() const;
  void setVariance(double val)
  {
    variance = val;
  }
  void setInitParam();
  double diagComputeElement(const CMatrix& X, unsigned int index) const;
  void setParam(double val, unsigned int paramNum);
  double getParam(unsigned int paramNum) const;
  void getGradX(CMatrix& g, const CMatrix& X, unsigned int pointNo, const CMatrix& X2, bool addG=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, bool addGrad=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, unsigned int index1, 
			const CMatrix& X2, unsigned int index2) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd, bool regularise=true) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd, bool regularise=true) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& X2, const CMatrix& cvGrd) const;
  double getGradParam(unsigned int index, const CMatrix& X, const CMatrix& cvGrd) const;
  void writeParamsToStream(ostream& out) const;
  void readParamsFromStream(istream& in);
  
  
 private:
  void _init();
  double weightVariance;
  double biasVariance;
  double variance;
  double degree;
  mutable CMatrix innerProdi;
  mutable CMatrix innerProdj;
  mutable CMatrix gscales;
};

ostream& operator<<(ostream& os, const CKern& A);
void writeKernToStream(const CKern& kern, ostream& out);
CKern* readKernFromStream(istream& in);




#endif
