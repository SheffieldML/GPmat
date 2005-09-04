#ifndef CDATAMODEL_H
#define CDATAMODEL_H
#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "ndlexceptions.h"
#include "ndlstrutil.h"
#include "COptimisable.h"
#include "CMatlab.h"
#include "CMatrix.h"
#include "CKern.h"
#include "CNoise.h"

class CDataModel : public CMatinterface // a model of a data set.
{
 public:
  // Initialise the model.
  virtual void init()
    {
      initStoreage();
      initVals();
    }
  // Initialise the storeage for the model.
  virtual void initStoreage()=0;
  // Set the initial values for the model.
  virtual void initVals()=0;
  virtual void display(ostream& os) const=0;
#ifdef _NDLMATLAB
  virtual mxArray* toMxArray() const=0;
  virtual void fromMxArray(const mxArray* matlabArray)=0;
#endif

};

/// This should form the base class for kernel models.
class CKernelModel : public CDataModel
{
 public:
  // The number of processes associated with the kernel.
  virtual int getNumProcesses() const=0;
  // the number of points used for computing core kernel.
  virtual int getNumActiveData() const=0;
  // The value of the inverse `jitter' for that point
  virtual double getBetaVal(const int i, const int j) const=0;

  // Update the stored kernel matrix.
  void updateK() const
    {
      double kVal=0.0;
      for(int i=0; i<getNumActiveData(); i++)
	{
	  K.setVal(kern.diagComputeElement(activeX, i), i, i);
	  for(int j=0; j<i; j++)
	    {
	      kVal=kern.computeElement(activeX, i, activeX, j);
	      K.setVal(kVal, i, j);
	      K.setVal(kVal, j, i);
	    }
	}
      K.setSymmetric(true);
    }
  
  // Update the stored inverse kernel matrix.
  void updateInvK(const int dim) const
    {
      invK.deepCopy(K);
      for(int i=0; i<getNumActiveData(); i++)
	invK.setVal(invK.getVal(i, i) + 1/getBetaVal(i, dim), i, i);
      invK.setSymmetric(true);
      CMatrix U(chol(invK));
      logDetK = invK.logDet(U); 
      invK.pdinv(U);
      
    }

  CKern& kern;
  CMatrix activeX;

  mutable CMatrix K;
  mutable CMatrix invK;
  mutable double logDetK;
};
class COptimisableModel : public CDataModel, public COptimisable
{
 public:
  virtual inline void setVerbosity(const int val) const
    {
      verbosity = val;
    }  
  virtual inline int getVerbosity() const
    {
      return verbosity;
    }
 private:
  mutable int verbosity;
  
};
class COptimisableKernelModel : public CKernelModel, public COptimisable
{
 public:
  virtual inline void setVerbosity(const int val) const
    {
      verbosity = val;
    }  
  virtual inline int getVerbosity() const
    {
      return verbosity;
    }
 private:
  mutable int verbosity;
  
};
class CMapModel : public CDataModel // a model which maps from one data space to another.
{
 public:
  // Initialise the storeage for the model.
  virtual void initStoreage()=0;
  // Set the initial values for the model.
  virtual void initVals()=0;
  virtual void display(ostream& os) const=0;
  virtual void toStream(ostream& out)=0;
  virtual void fromStream(const istream& out) const=0;
#ifdef _NDLMATLAB
  virtual mxArray* toMxArray() const=0;
  virtual void fromMxArray(const mxArray* matlabArray)=0;
#endif
 private:
};

#endif
