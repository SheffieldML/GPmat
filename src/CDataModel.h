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
