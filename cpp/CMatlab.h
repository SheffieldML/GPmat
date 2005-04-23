#ifndef CMATLAB_H
#define CMATLAB_H

#include "mex.h"
#include "mat.h"
using namespace std;

class CMatinterface {
 public:
  virtual mxArray* toMxArray() const=0;
  virtual void fromMxArray(const mxArray* matlabArray)=0;
  void readMatlabFile(const string fileName, const string variableName) 
    {
      MATFile* matFile = matOpen(fileName.c_str(), "r");
      if(matFile==NULL)
	{
	  cerr << "Unable to open file " << fileName << endl;
	}
      
      mxArray* matlabArray = matGetVariable(matFile, variableName.c_str());
      if(matlabArray==NULL)
	{
	  cerr << "Variable " << variableName << " not present in " << fileName << endl;
	}
      fromMxArray(matlabArray);
      mxDestroyArray(matlabArray);
      if(matClose(matFile) !=0 )
	{
	  cerr << "Unable to close file " << fileName << endl;
	}
      
    }
  void updateMatlabFile(string fileName, const string variableName) const
    {
      MATFile* matFile = matOpen(fileName.c_str(), "u");
      if(matFile==NULL)
	{      
	  cerr << "Unable to open file " << fileName << endl;
	}
      mxArray* matlabArray = toMxArray();
      matPutVariable(matFile, variableName.c_str(), matlabArray);
      mxDestroyArray(matlabArray);
      if(matClose(matFile) !=0 )
	{
	  cerr << "Unable to close file " << fileName << endl;
	}
    }
    
  void writeMatlabFile(const string fileName, const string variableName) const
    {
      MATFile* matFile = matOpen(fileName.c_str(), "w");
      if(matFile==NULL)
	{      
	  cerr << "Unable to open file " << fileName << endl;
	}
      mxArray* matlabArray = toMxArray();
      matPutVariable(matFile, variableName.c_str(), matlabArray);
      mxDestroyArray(matlabArray);
      if(matClose(matFile) !=0 )
	{
	  cerr << "Unable to close file " << fileName << endl;
	}
    }
  mxArray* convertMxArray(const bool val) const
    {
      int dims[1];
      dims[0] = 1;
      mxArray* matlabArray = mxCreateNumericArray(1, dims, mxDOUBLE_CLASS, mxREAL);
      double* matlabVals = mxGetPr(matlabArray);
      if(val)
	matlabVals[0] = 1.0;
      else
	matlabVals[0] = 0.0;
      return matlabArray;
    }
  mxArray* convertMxArray(const double val) const
    {
      int dims[1];
      dims[0] = 1;
      mxArray* matlabArray = mxCreateNumericArray(1, dims, mxDOUBLE_CLASS, mxREAL);
      double* matlabVals = mxGetPr(matlabArray);
      matlabVals[0] = val;
      return matlabArray;
    }
  mxArray* convertMxArray(const string val) const
    {
      return mxCreateString(val.c_str());
    }
  mxArray* convertMxArray(const vector<int> vals) const
    {
      int dims[2];
      dims[0]=vals.size();
      dims[1] = 1;
      mxArray* matlabArray = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      double* matlabVals = mxGetPr(matlabArray);
      for(int i=0; i<vals.size(); i++)
	{
	  matlabVals[i]=(double)vals[i];
	}
      return matlabArray;
    }
  int mxArrayToInt(const mxArray* matlabArray) const
    {
      int val = 0;
      mxClassID classID = mxGetClassID(matlabArray);
      if(classID==mxDOUBLE_CLASS)
	{
	  assert(mxGetN(matlabArray)==1);
	  assert(mxGetM(matlabArray)==1);
	  double* valD = mxGetPr(matlabArray);
	  val = (int)valD[0];
	}
      else
	cerr << "Conversion of this type to int not yet supported." << endl;
      return val;
    }
  double mxArrayToDouble(const mxArray* matlabArray) const
    {
      double val = 0;
      mxClassID classID = mxGetClassID(matlabArray);
      if(classID==mxDOUBLE_CLASS)
	{
	  assert(mxGetN(matlabArray)==1);
	  assert(mxGetM(matlabArray)==1);
	  double* valD = mxGetPr(matlabArray);
	  val = valD[0];
	}
      else
	  cerr << "Conversion of this type to double not yet supported." << endl;
      return val;
    }
  string mxArrayToString(const mxArray* matlabArray) const
    {
      char* charVal;
      int buflen;
      mxClassID classID = mxGetClassID(matlabArray);
      if(classID ==  mxCHAR_CLASS)
	{
	  buflen=mxGetNumberOfElements(matlabArray)*sizeof(char) + 1;
	  int ret = mxGetString(matlabArray, charVal, buflen);
	}
      else
	cerr << "Conversion of this type to string is not yet supported." << endl;
      string val(charVal);
      return val;
    }
  vector<int> mxArrayToVectorInt(const mxArray* matlabArray) const
    {
      vector<int> val;
      mxClassID classID = mxGetClassID(matlabArray);
      if(classID == mxDOUBLE_CLASS)
	{
	  int length=mxGetNumberOfElements(matlabArray);
	  double* valD = mxGetPr(matlabArray);
	  for(int i=0; i<length; i++)
	    val.push_back((int)valD[i]);
	}
      else
	cerr << "Conversion of this type to vector<int> is not yet supported." << endl;
      return val;
    }
  bool mxArrayToBool(const mxArray* matlabArray) const
    {
      bool val;
      mxClassID classID = mxGetClassID(matlabArray);
      if(classID==mxDOUBLE_CLASS)
	{
	  assert(mxGetN(matlabArray)==1);
	  assert(mxGetM(matlabArray)==1);
	  double* valD = mxGetPr(matlabArray);
	  val = valD[0];
	  if(val!=0)
	    val=true;
	  else
	    val=false;
	}
      else
	cerr << "Conversion of this type to double not yet supported." << endl;
    }	  
  
  
  mxArray* mxArrayExtractMxArrayField(const mxArray* matlabArray, const string fieldName, const int index=0) const
    {
      const char* fName = fieldName.c_str();
      if(mxGetClassID(matlabArray) != mxSTRUCT_CLASS)
	cerr << "mxArray is not a structure." << endl;
      mxArray* fieldPtr = mxGetField(matlabArray, index, fName);
      if(fieldPtr == NULL)
	cout << "No such field as " << fieldName << "." << endl;
      return fieldPtr;
    }
	
  int mxArrayExtractIntField(const mxArray* matlabArray, const string fieldName, const int index=0) const
    {
      mxArray* fieldPtr = mxArrayExtractMxArrayField(matlabArray, fieldName, index);
      return mxArrayToInt(fieldPtr);
    }
  double mxArrayExtractDoubleField(const mxArray* matlabArray, const string fieldName, const int index=0) const
    {
      mxArray* fieldPtr = mxArrayExtractMxArrayField(matlabArray, fieldName, index);
      return mxArrayToDouble(fieldPtr);
    }
  string mxArrayExtractStringField(const mxArray* matlabArray, const string fieldName, const int index=0) const
    {
      mxArray* fieldPtr = mxArrayExtractMxArrayField(matlabArray, fieldName, index);
      return mxArrayToString(fieldPtr);
    }

  vector<int> mxArrayExtractVectorIntField(const mxArray* matlabArray, const string fieldName, const int index=0) const
    {
      mxArray* fieldPtr = mxArrayExtractMxArrayField(matlabArray, fieldName, index);
      return mxArrayToVectorInt(fieldPtr);
    }
  bool mxArrayExtractBoolField(const mxArray* matlabArray, const string fieldName, const int index=0)
    {
      mxArray* fieldPtr = mxArrayExtractMxArrayField(matlabArray, fieldName, index);
      return mxArrayToBool(fieldPtr);
    }

};

#endif
