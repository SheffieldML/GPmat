#ifndef CKERN_H
#define CKERN_H
#include <cmath>
#include "CTransform.h"
#include "CMatrix.h"
#include "CDist.h"
#include "ndlstrutil.h"
#include "ndlutil.h"
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <algorithm>
using namespace std;

const string KERNVERSION="0.1";

// The base class for all kernels. It implements: CMatinterface, for
// saving and loading from MATLAB; CRegularisable, for placing priors
// over parameters; CTransformable, for allowing parameters to be
// transformed so that they are only optimised in, for example,
// positive half-spaces
class CKern : public CMatinterface, public CTransformable, public CRegularisable {
 public:
  CKern()
    {}
  CKern(const CMatrix& X){}
  
  virtual ~CKern(){}
  CKern(const CKern& kern){}
  virtual CKern* clone() const=0;
  // set initial parameters.
  virtual void setInitParam()=0;
  // compute an element of the diagonal.
  virtual double diagComputeElement(const CMatrix& X, int index) const=0;
  // compute the entire diagonal
  virtual void diagCompute(CMatrix& d, const CMatrix& X) const
    {
      assert(X.rowsMatch(d));
      assert(d.getCols()==1);
      for(int i=0; i<X.getRows(); i++)
	d.setVal(diagComputeElement(X, i), i);
    }
  // Compute the diagonal at particualr indices.
  virtual void diagCompute(CMatrix& d, const CMatrix& X, const vector<int> indices) const
    {
      assert(d.getRows()==indices.size());
      assert(d.getCols()==1);
      for(int i=0; i<indices.size(); i++)
	d.setVal(diagComputeElement(X, indices[i]), i);
    }
  // Set the parameters of the kernel.
  virtual void setParam(double, int)=0;
  // Get gradients of the kernel with respect to input values.
  //   g[i].val(k,j) = d kern(X_row_i,X2_row_k)/ d x_component_j
  virtual void getGradX(vector<CMatrix*>& g, const CMatrix& X, const CMatrix& X2, const bool addG=false) const=0;
  // Get gradients of the kernel diagonal with respect to input values.
  // WVB: I think this is separate from getGradX just because the diagonal vals
  // (i.e. kern(X_row_i,X_row_i) should always be taken from the same matrix X.
  // Not to imply that I understand why we need X and X2 in getGradX, really.
  virtual void getDiagGradX(CMatrix& g, const CMatrix& X, const bool addG=false) const=0;
  // Return the white noise component of the kernel.
  virtual double getWhite() const
    {
      return 0.0;
    }
  // Compute an element of the kernel matrix.
  virtual double computeElement(const CMatrix& X1, int index1,
			 const CMatrix& X2, int index2) const=0;
  // Compute specified rows and columns of the kernel matrix.
  virtual void compute(CMatrix& K, const CMatrix& X1, const vector<int> indices1,
			  const CMatrix& X2, const vector<int> indices2) const
    {
      assert(K.rowsMatch(X1));
      assert(K.getCols()==X2.getRows());
      for(int i=0; i<indices1.size(); i++)
	{
	  for(int j=0; j<indices2.size(); j++)
	    {
	      assert(indices1[i]<K.getRows());
	      assert(indices2[j]<K.getCols());
	      K.setVal(computeElement(X1, indices1[i], X2, indices2[j]), indices1[i], indices2[j]);
	    }
	}
    }
  // Compute the kernel matrix for a data set.
  virtual void compute(CMatrix& K, const CMatrix& X) const
    {
      assert(K.rowsMatch(X));
      assert(K.isSquare());
      double k = 0.0;
      for(int i=0; i<K.getRows(); i++)
	{
	  for(int j=0; j<i; j++)
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
      assert(K.rowsMatch(X));
      assert(K.getCols()==X2.getRows());
      for(int i=0; i<K.getRows(); i++)
	{
	  for(int j=0; j<K.getCols(); j++)
	    {
	      K.setVal(computeElement(X, i, X2, j), i, j);
	    }
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
  virtual void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd) const
    {
      assert(g.getRows()==1);
      assert(g.getCols()==nParams);
      assert(X.rowsMatch(cvGrd));
      assert(cvGrd.isSquare());
      for(int i=0; i<nParams; i++)
	g.setVal(getGradParam(i, X, cvGrd), i);
      addPriorGrad(g); /// don't forget to add prior gradient at the end.
    }
  // Get gradient of a particular parameter.
  virtual double getGradParam(int index, const CMatrix& X, const CMatrix& cvGrd) const=0;
  // Get a particular parameter.
  // For compound kernels when a new kernel is added.
  virtual int addKern(CKern* kern)
    {
      cerr << "You cannot add a kernel to this kernel." << endl;
      return 0;
    }
  virtual double getParam(int) const=0;
  // Set the parameters from a vector of parameters.
  void setParams(const CMatrix& paramVec)
    {
      for(int i=0; i<nParams; i++)
	setParam(paramVec.getVal(i), i);
    }
  // Place the parameters in a vector.
  void getParams(CMatrix& paramVec) const
    {
      for(int i=0; i<nParams; i++)
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
  inline void setInputDim(int dim)
    {
      inputDim = dim;
      setInitParam();
    }
  // Get the input dimension.
  inline int getInputDim() const
    {
      return inputDim;
    }
  // How many kernel parameters are there?
  inline int getNumParams() const
    {
      return nParams;
    }
  // Assign a name to the kernel parameters.
  void setParamName(const string name, int index)
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
  // Get the name of a kernel parameter.
  virtual string getParamName(int index) const
    {
      assert(index>=0);
      assert(index<paramNames.size());
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
  ostream& display(ostream& os) const;
  
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
  void getGradTransParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd) const;
  // specify tests for equality between kernels.
  bool equals(const CKern& kern, double tol=ndlutil::MATCHTOL) const;
  
 protected:
  int nParams;
  string kernName;
  string type;
  vector<string> paramNames;
 private:
  int inputDim;
};

// CArdKern is the base class for any kernel that uses multiple input parameters.
class CArdKern : public CKern {
#ifdef _NDLMATLAB
  virtual void addParamToMxArray(mxArray* matlabArray) const;
  // Gets the parameters from the mxArray.
  virtual void extractParamFromMxArray(const mxArray* matlabArray);
  // returns sum(sum(cvGrd.*dK/dparam)) 
#endif
 protected: 
  CMatrix scales;
};

// Compound Kernel --- This kernel combines other kernels additively together.
class CCmpndKern: public CKern {
 public:
  CCmpndKern();
  CCmpndKern(int inDim);
  CCmpndKern(const CMatrix& X);
  ~CCmpndKern();
  CCmpndKern(const CCmpndKern&);
  CCmpndKern* clone() const
    {
      return new CCmpndKern(*this);
    }
  CCmpndKern(vector<CKern*> kernels);
  
  void setInitParam();
  double diagComputeElement(const CMatrix& X, int index1) const;
  void diagCompute(CMatrix& d, const CMatrix& X) const;
  void setParam(double val, int paramNum);
  double getParam(int paramNum) const;
  string getParamName(int index) const;
  void getGradX(vector<CMatrix*>& g, const CMatrix& X, const CMatrix& X2, const bool addG=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, const bool addG=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, int index1, 
		 const CMatrix& X2, int index2) const;
  void compute(CMatrix& K, const CMatrix& X, const CMatrix& X2) const;
  void compute(CMatrix& K, const CMatrix& X) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd) const;
  double getGradParam(int index, const CMatrix& X, const CMatrix& cvGrd) const;
  double priorLogProb() const;
  void addPrior(CDist* prior, int index);


  void writeParamsToStream(ostream& out) const;
  void readParamsFromStream(istream& in);
  int addKern(CKern* kern);

  int getNumKerns() const
    {
      return components.size();
    }
#ifdef _NDLMATLAB
  // Kernel specific code for toMxArray() call.
  void addParamToMxArray(mxArray* matlabArray) const;
  // Kernel specific code for fromMxArray() call.
  void extractParamFromMxArray(const mxArray* matlabArray);
#endif
 private:
  // this is a heterogeneous container.
  vector<CKern*> components;
  

};
// White Noise Kernel.
class CWhiteKern: public CKern {
 public:
  CWhiteKern();
  CWhiteKern(int inDim);
  CWhiteKern(const CMatrix& X);
  ~CWhiteKern();
  CWhiteKern(const CWhiteKern&);
  CWhiteKern* clone() const
    {
      return new CWhiteKern(*this);
    }
 void setInitParam();
  double diagComputeElement(const CMatrix& X, int index) const;
  void diagCompute(CMatrix& d, const CMatrix& X) const;
  void setParam(double val, int paramNum);
  double getParam(int paramNum) const;
  void getGradX(vector<CMatrix*>& g, const CMatrix& X, const CMatrix& X2, const bool addG=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, const bool addG=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, int index1, 
		 const CMatrix& X2, int index2) const;
  void  compute(CMatrix& K, const CMatrix& X, const CMatrix& X2) const;
  void compute(CMatrix& K, const CMatrix& X) const;
  double getGradParam(int index, const CMatrix& X, const CMatrix& cvGrd) const;


 private:
  double variance;
};

// Bias Kernel.
class CBiasKern: public CKern {
 public:
  CBiasKern();
  CBiasKern(int inDim);
  CBiasKern(const CMatrix& X);
  ~CBiasKern();
  CBiasKern(const CBiasKern&);
  CBiasKern* clone() const
    {
      return new CBiasKern(*this);
    }
  void setInitParam();
  double diagComputeElement(const CMatrix& X, int index) const;
  void diagCompute(CMatrix& d, const CMatrix& X) const;
  void setParam(double val, int paramNum);
  double getParam(int paramNum) const;
  void getGradX(vector<CMatrix*>& g, const CMatrix& X, const CMatrix& X2, const bool addG=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, const bool addG=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, int index1,
		 const CMatrix& X2, int index2) const;
  void compute(CMatrix& K, const CMatrix& X, const CMatrix& X2) const;
  void compute(CMatrix& K, const CMatrix& X) const;
  double getGradParam(int index, const CMatrix& X, const CMatrix& cvGrd) const;



 private:
  double variance;

};
// RBF Kernel, also known as the Gaussian or squared exponential kernel.
class CRbfKern: public CKern {
 public:
  CRbfKern();
  CRbfKern(int inDim);
  CRbfKern(const CMatrix& X);
  ~CRbfKern();
  CRbfKern(const CRbfKern&);
  CRbfKern* clone() const
    {
      return new CRbfKern(*this);
    }
  void setInitParam();
  double diagComputeElement(const CMatrix& X, int index) const;
  void diagCompute(CMatrix& d, const CMatrix& X) const;
  void setParam(double val, int paramNum);
  double getParam(int paramNum) const;
  void getGradX(vector<CMatrix*>& g, const CMatrix& X, const CMatrix& X2, const bool addG=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, const bool addG=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, int index1, 
		  const CMatrix& X2, int index2) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd) const;
  double getGradParam(int index, const CMatrix& X, const CMatrix& cvGrd) const;

 private:
  double variance;
  double inverseWidth;

};

// Linear Kernel, also known as the inner product kernel.
class CLinKern: public CKern {
 public:
  CLinKern();
  CLinKern(int inDim);
  CLinKern(const CMatrix& X);
  ~CLinKern();
  CLinKern(const CLinKern&);
  CLinKern* clone() const
    {
      return new CLinKern(*this);
    }
  void setInitParam();
  double diagComputeElement(const CMatrix& X, int index) const;
  void setParam(double val, int paramNum);
  double getParam(int paramNum) const;
  void getGradX(vector<CMatrix*>& g, const CMatrix& X, const CMatrix& X2, const bool addG=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, const bool addG=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, int index1, 
		 const CMatrix& X2, int index2) const;
  void compute(CMatrix& K, const CMatrix& X, const CMatrix& X2) const;
  void compute(CMatrix& K, const CMatrix& X) const;
  double getGradParam(int index, const CMatrix& X, const CMatrix& cvGrd) const;


 private:
  double variance;
};

// MLP Kernel or arcsin kernel. Based on a multi-layer perceptron with infinite hidden nodes. See Williams (1996) "Computing with Infinite Networks" in NIPS 9.
class CMlpKern: public CKern {
 public:
  CMlpKern();
  CMlpKern(int inDim);
  CMlpKern(const CMatrix& X);
  ~CMlpKern();
  CMlpKern(const CMlpKern&);
  CMlpKern* clone() const
    {
      return new CMlpKern(*this);
    }
  void setInitParam();
  double diagComputeElement(const CMatrix& X, int index) const;
  void setParam(double val, int paramNum);
  double getParam(int paramNum) const;
  void getGradX(vector<CMatrix*>& g, const CMatrix& X, const CMatrix& X2, const bool addG=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, const bool addG=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, int index1, 
		  const CMatrix& X2, int index2) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd) const;
  double getGradParam(int index, const CMatrix& X, const CMatrix& cvGrd) const;

 private:
  double weightVariance;
  double biasVariance;
  double variance;
  mutable CMatrix innerProd;
};

// Polynomial Kernel, not generally recommended as it `extreme behaviour' outside the region where the argument's absolute value is less than 1.
class CPolyKern: public CKern {
 public:
  CPolyKern();
  CPolyKern(int inDim);
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
  void setInitParam();
  double diagComputeElement(const CMatrix& X, int index) const;
  void setParam(double val, int paramNum);
  double getParam(int paramNum) const;
  void getGradX(vector<CMatrix*>& g, const CMatrix& X, const CMatrix& X2, const bool addG=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, const bool addG=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, int index1, 
		  const CMatrix& X2, int index2) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd) const;
  double getGradParam(int index, const CMatrix& X, const CMatrix& cvGrd) const;
  void writeParamsToStream(ostream& out) const;
  void readParamsFromStream(istream& in);
#ifdef _NDLMATLAB
  void addParamToMxArray(mxArray* matlabArray) const;
  void extractParamFromMxArray(const mxArray* matlabArray);
#endif
 private:
  double weightVariance;
  double biasVariance;
  double variance;
  double degree;
  mutable CMatrix innerProd;
};

// Linear ARD Kernel --- automatic relevance determination version of the linear kernel.
class CLinardKern: public CArdKern {
 public:
  CLinardKern();
  CLinardKern(int inDim);
  CLinardKern(const CMatrix& X);
  ~CLinardKern();
  CLinardKern(const CLinardKern&);
  CLinardKern* clone() const
    {
      return new CLinardKern(*this);
    }
  void setInitParam();
  double diagComputeElement(const CMatrix& X, int index) const;
  void setParam(double val, int paramNum);
  double getParam(int paramNum) const;
  void getGradX(vector<CMatrix*>& g, const CMatrix& X, const CMatrix& X2, const bool addGrad=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, const bool addGrad=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, int index1, 
		 const CMatrix& X2, int index2) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd) const;
  double getGradParam(int index, const CMatrix& X, const CMatrix& cvGrd) const;


 private:
  double variance;
};

// RBF ARD Kernel --- automatic relevance determination of the RBF kernel.
class CRbfardKern: public CArdKern {
 public:
  CRbfardKern();
  CRbfardKern(int inDim);
  CRbfardKern(const CMatrix& X);
  ~CRbfardKern();
  CRbfardKern(const CRbfardKern&);
  CRbfardKern* clone() const
    {
      return new CRbfardKern(*this);
    }
  void setInitParam();
  double diagComputeElement(const CMatrix& X, int index) const;
  void setParam(double val, int paramNum);
  double getParam(int paramNum) const;
  void getGradX(vector<CMatrix*>& g, const CMatrix& X, const CMatrix& X2, const bool addGrad=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, const bool addGrad=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, int index1, 
		 const CMatrix& X2, int index2) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd) const;
  double getGradParam(int index, const CMatrix& X, const CMatrix& cvGrd) const;


 private:
  double variance;
  double inverseWidth;
  mutable CMatrix gscales;
};

// MLP ARD Kernel --- automatic relevance determination version of the MLP kernel.
class CMlpardKern: public CArdKern {
 public:
  CMlpardKern();
  CMlpardKern(int inDim);
  CMlpardKern(const CMatrix& X);
  ~CMlpardKern();
  CMlpardKern(const CMlpardKern&);
  CMlpardKern* clone() const
    {
      return new CMlpardKern(*this);
    }
  void setInitParam();
  double diagComputeElement(const CMatrix& X, int index) const;
  void setParam(double val, int paramNum);
  double getParam(int paramNum) const;
  void getGradX(vector<CMatrix*>& g, const CMatrix& X, const CMatrix& X2, const bool addGrad=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, const bool addGrad=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, int index1, 
		 const CMatrix& X2, int index2) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd) const;
  double getGradParam(int index, const CMatrix& X, const CMatrix& cvGrd) const;


 private:
  double weightVariance;
  double biasVariance;
  double variance;
  mutable CMatrix innerProd;
  mutable CMatrix gscales;
};

// Polynomial ARD Kernel --- automatic relevance determination version of the polynomial kernel.
class CPolyardKern: public CArdKern {
 public:
  CPolyardKern();
  CPolyardKern(int inDim);
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
  void setInitParam();
  double diagComputeElement(const CMatrix& X, int index) const;
  void setParam(double val, int paramNum);
  double getParam(int paramNum) const;
  void getGradX(vector<CMatrix*>& g, const CMatrix& X, const CMatrix& X2, const bool addGrad=false) const;
  void getDiagGradX(CMatrix& g, const CMatrix& X, const bool addGrad=false) const;
  double getWhite() const;
  double computeElement(const CMatrix& X1, int index1, 
		 const CMatrix& X2, int index2) const;
  void getGradParams(CMatrix& g, const CMatrix& X, const CMatrix& cvGrd) const;
  double getGradParam(int index, const CMatrix& X, const CMatrix& cvGrd) const;
  void writeParamsToStream(ostream& out) const;
  void readParamsFromStream(istream& in);


 private:
  double weightVariance;
  double biasVariance;
  double variance;
  double degree;
  mutable CMatrix innerProd;
  mutable CMatrix gscales;
};

ostream& operator<<(ostream& os, const CKern& A);
void writeKernToStream(const CKern& kern, ostream& out);
CKern* readKernFromStream(istream& in);


#endif
