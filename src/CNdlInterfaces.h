#ifndef CNdlInterfaces_H
#define CNdlInterfaces_H

#include <string>
#include <iostream>
#include <sstream>

#include "ndlexceptions.h"
#include "ndlstrutil.h"

#ifdef _NDLMATLAB
#include "mex.h"
#include "mat.h"
#endif

const double MINVERSION=0.2;
const double VERSION=0.2;

using namespace std;

class CStreamInterface 
{
public:
  virtual ~CStreamInterface() {};
  virtual void toStream(ostream& out) const
  {
    out << setiosflags(ios::fixed);
    out << setprecision(6);
    writeToStream(out, "version", getCurrentVersion());
    out << setiosflags(ios::scientific);
    out << setprecision(17);
    writeParamsToStream(out);
  }
  static double readVersionFromStream(istream& in) throw(ndlexceptions::StreamVersionError&)
  {
    double ver = readDoubleFromStream(in, "version");
    if(ver<getMinCompatVersion())
      throw ndlexceptions::StreamVersionError();
    return ver;

  }
  static int readIntFromStream(istream& in, const std::string fieldName)
  {
    string str = readStringFromStream(in, fieldName);
    return atol(str.c_str());
  }
  static double readDoubleFromStream(istream& in, const std::string fieldName)
  {
    string str = readStringFromStream(in, fieldName);
    return atof(str.c_str());
}
  static bool readBoolFromStream(istream& in, const std::string fieldName)
  {
    string str = readStringFromStream(in, fieldName);
    if(atol(str.c_str())!=0)
      return true;
    else
      return false;
  }
  static vector<int> readVectorIntFromStream(istream& in, const std::string fieldName) throw(ndlexceptions::StreamFormatError&)
  {
    vector<int> vals;
    string str = readStringFromStream(in, fieldName);
    vector<string> tokens;
    ndlstrutil::tokenise(tokens, str, " ");
    if(tokens.size()==0)
      throw ndlexceptions::StreamFormatError(fieldName, "Zero length vector<int>.");
    for(size_t i=0; i<tokens.size(); i++)
      vals.push_back(atol(tokens[i].c_str()));
    return vals;
  }
  static vector<unsigned int> readVectorUintFromStream(istream& in, const std::string fieldName) throw(ndlexceptions::StreamFormatError&)
  {
    vector<unsigned int> vals;
    string str = readStringFromStream(in, fieldName);
    vector<string> tokens;
    ndlstrutil::tokenise(tokens, str, " ");
    if(tokens.size()==0)
      throw ndlexceptions::StreamFormatError(fieldName, "Zero length vector<unsigned int>.");
    for(size_t i=0; i<tokens.size(); i++)
      vals.push_back(atol(tokens[i].c_str()));
    return vals;
  }
  static string readStringFromStream(istream& in, const std::string fieldName) throw(ndlexceptions::StreamFormatError&)
  {
    string line;
    vector<string> tokens;
    ndlstrutil::getline(in, line);
    ndlstrutil::tokenise(tokens, line, "=");
    if(tokens.size()!=2 || tokens[0]!=fieldName)
      throw ndlexceptions::StreamFormatError(fieldName);
    return tokens[1];
}
  static void writeToStream(ostream& out, const std::string fieldName, const vector<int> val)
  {
    out << fieldName << "=";
    for(size_t i=0; i<val.size()-1; i++)
      out << ndlstrutil::itoa(val[i]) << " ";
    out << ndlstrutil::itoa(val[val.size()-1]) << endl;
  }
  static void writeToStream(ostream& out, const std::string fieldName, const vector<unsigned int> val)
  {
    out << fieldName << "=";
    for(size_t i=0; i<val.size()-1; i++)
      out << ndlstrutil::itoa(val[i]) << " ";
    out << ndlstrutil::itoa(val[val.size()-1]) << endl;
  }
  static void writeToStream(ostream& out, const std::string fieldName, const int val)
  {
    out << fieldName << "=" << ndlstrutil::itoa(val) << endl; 
  }
  static void writeToStream(ostream& out, const std::string fieldName, const unsigned int val)
  {
    out << fieldName << "=" << ndlstrutil::itoa(val) << endl; 
  }
  static void writeToStream(ostream& out, const std::string fieldName, const double val)
  {
    out << fieldName << "=" << val << endl;
  }
  static void writeToStream(ostream& out, const std::string fieldName, const bool val)
  {
    out << fieldName << "=";
    if(val)
      out << "1" << endl;
    else
      out << "0" << endl;
  }
  static void writeToStream(ostream& out, const std::string fieldName, const std::string val)
  {
      out << fieldName << "=" << val << endl;
  }
  static void writeToStream(ostream& out, const std::string fieldName, const char* val)
  {
      out << fieldName << "=" << val << endl;
  }
  
  static string getBaseTypeStream(istream& in)
  {
    return readStringFromStream(in, "baseType");
  }
  static string getTypeStream(istream& in)
  {
    return readStringFromStream(in, "type");
  }
  virtual void fromStream(istream& in) throw(ndlexceptions::StreamFormatError&)
  {
    readVersionFromStream(in);
    readParamsFromStream(in);
  }

  double getCurrentVersion() const 
  {
    return VERSION;
  }
  static double getMinCompatVersion()
  {
    return MINVERSION;
  }
  virtual void writeParamsToStream(ostream& out) const=0;
  virtual void readParamsFromStream(istream& out)=0;
  void toFile(const string fileName, const string comment="") const throw(ndlexceptions::FileWriteError&)
  {
    ofstream out(fileName.c_str());
    if(!out) throw ndlexceptions::FileWriteError(fileName);
    if(comment.size()>0)
      out << "# " << comment << endl;
    toStream(out);
    out.close();
  }
  void fromFile(const string fileName) throw(ndlexceptions::FileReadError&, ndlexceptions::FileFormatError&)
  {
    ifstream in(fileName.c_str());
    if(!in.is_open()) throw ndlexceptions::FileReadError(fileName);
    try
    {
      fromStream(in);
    }
    catch(ndlexceptions::StreamFormatError& err)
    {
      throw ndlexceptions::FileFormatError(fileName, err);
    } 
    in.close();
  }
    
    #ifdef _HDF5
        virtual void writeToHdf5( const std::string& filename, const std::string& path_to_dataset ) const=0;
        virtual void readFromHdf5( const std::string& filename, const std::string& path_to_dataset )=0;
    #endif
};

#ifdef _NDLMATLAB

// An abstract base class which enables loading and saving of a class to MATLAB.
class CMatInterface 
{
 public:
  virtual mxArray* toMxArray() const=0;
  virtual void fromMxArray(const mxArray* matlabArray)=0;
  void readMatlabFile(const string fileName, const string variableName) throw(ndlexceptions::FileReadError&, ndlexceptions::FileFormatError&)
  {
    MATFile* matFile = matOpen(fileName.c_str(), "r");
    if(matFile==NULL)
      throw ndlexceptions::FileReadError(fileName);
    mxArray* matlabArray = matGetVariable(matFile, variableName.c_str());
    if(matlabArray==NULL)
    {
      ndlexceptions::MatlabInterfaceReadError matError(variableName);
      throw ndlexceptions::FileFormatError(fileName, matError);
    }
    try
    {
      fromMxArray(matlabArray);
    }
    catch(ndlexceptions::MatlabInterfaceReadError& err)
    {
      throw ndlexceptions::FileFormatError(fileName, err);
    }
    mxDestroyArray(matlabArray);
    if(matClose(matFile) !=0 )
      throw ndlexceptions::FileReadError(fileName);   
  }
  void updateMatlabFile(string fileName, const string variableName) const throw(ndlexceptions::FileWriteError&)
  {
    MATFile* matFile = matOpen(fileName.c_str(), "u");
    if(matFile==NULL)
      throw ndlexceptions::FileWriteError(fileName);
    mxArray* matlabArray = toMxArray();
    matPutVariable(matFile, variableName.c_str(), matlabArray);
    mxDestroyArray(matlabArray);
    if(matClose(matFile) !=0 )
      throw ndlexceptions::FileWriteError(fileName);
  }
  
  void writeMatlabFile(const string fileName, const string variableName) const throw(ndlexceptions::FileWriteError&)
  {
    MATFile* matFile = matOpen(fileName.c_str(), "w");
    if(matFile==NULL)
      throw ndlexceptions::FileWriteError(fileName);
    mxArray* matlabArray = toMxArray();
    matPutVariable(matFile, variableName.c_str(), matlabArray);
    mxDestroyArray(matlabArray);
    if(matClose(matFile) !=0 )
      throw ndlexceptions::FileWriteError(fileName);
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
    dims[0]=(int)vals.size();
    dims[1] = 1;
    mxArray* matlabArray = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    double* matlabVals = mxGetPr(matlabArray);
    for(size_t i=0; i<vals.size(); i++)
    {
      matlabVals[i]=(double)vals[i];
    }
      return matlabArray;
  }
  mxArray* convertMxArray(const vector<unsigned int> vals) const
  {
    int dims[2];
    dims[0]=(int)vals.size();
    dims[1] = 1;
    mxArray* matlabArray = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    double* matlabVals = mxGetPr(matlabArray);
    for(size_t i=0; i<vals.size(); i++)
    {
      matlabVals[i]=(double)vals[i];
    }
      return matlabArray;
  }
  int mxArrayToInt(const mxArray* matlabArray) const throw(ndlexceptions::NotImplementedError&)
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
      throw ndlexceptions::NotImplementedError("Conversion of this type to int not yet supported.");
    return val;
  }
  double mxArrayToDouble(const mxArray* matlabArray) const throw(ndlexceptions::NotImplementedError&)
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
      throw ndlexceptions::NotImplementedError("Conversion of this type to double not yet supported.");
    return val;
  }
  string mxArrayToString(const mxArray* matlabArray) const throw(ndlexceptions::NotImplementedError&)
  {
    vector<char> charVal;
    int buflen;
    mxClassID classID = mxGetClassID(matlabArray);
    if(classID ==  mxCHAR_CLASS)
    {
      buflen=mxGetNumberOfElements(matlabArray)*sizeof(char) + 1;
      charVal.resize(buflen);
      int ret = mxGetString(matlabArray, &charVal[0], buflen);
    }
    else
      throw ndlexceptions::NotImplementedError("Conversion of this type to string not yet supported.");
    string val(&charVal[0]);
    return val;
  }
  vector<int> mxArrayToVectorInt(const mxArray* matlabArray) const throw(ndlexceptions::NotImplementedError&)
  {
    vector<int> val;
    mxClassID classID = mxGetClassID(matlabArray);
    if(classID == mxDOUBLE_CLASS)
    {
      unsigned int length=mxGetNumberOfElements(matlabArray);
      double* valD = mxGetPr(matlabArray);
      for(unsigned int i=0; i<length; i++)
  val.push_back((int)valD[i]);
    }
    else
      throw ndlexceptions::NotImplementedError("Conversion of this type to vector<int> not yet supported.");
    return val;
  }
  vector<unsigned int> mxArrayToVectorUint(const mxArray* matlabArray) const throw(ndlexceptions::NotImplementedError&)
  {
    vector<unsigned int> val;
    mxClassID classID = mxGetClassID(matlabArray);
    if(classID == mxDOUBLE_CLASS)
    {
      unsigned int length=mxGetNumberOfElements(matlabArray);
      double* valD = mxGetPr(matlabArray);
      for(unsigned int i=0; i<length; i++)
  val.push_back((unsigned int)valD[i]);
    }
    else
      throw ndlexceptions::NotImplementedError("Conversion of this type to vector<int> not yet supported.");
    return val;
  }
  bool mxArrayToBool(const mxArray* matlabArray) const throw(ndlexceptions::NotImplementedError&)
  {
    bool val;
    mxClassID classID = mxGetClassID(matlabArray);
    if(classID==mxDOUBLE_CLASS)
    {
      assert(mxGetN(matlabArray)==1);
      assert(mxGetM(matlabArray)==1);
      double* valD = mxGetPr(matlabArray);
      if(valD!=0)
	val=true;
      else
	val=false;
    }
    else
      throw ndlexceptions::NotImplementedError("Conversion of this type to bool not yet supported.");
    return false;
  }	  
  
  
  mxArray* mxArrayExtractMxArrayField(const mxArray* matlabArray, const string fieldName, const int index=0) const throw(ndlexceptions::MatlabInterfaceReadError&)
  {
    const char* fName = fieldName.c_str();
    if(mxGetClassID(matlabArray) != mxSTRUCT_CLASS)
      throw ndlexceptions::MatlabInterfaceReadError(fieldName);
    mxArray* fieldPtr = mxGetField(matlabArray, index, fName);
    return fieldPtr;
  }
  
  int mxArrayExtractIntField(const mxArray* matlabArray, const string fieldName, const int index=0) const throw(ndlexceptions::MatlabInterfaceReadError&)
  {
    mxArray* fieldPtr = mxArrayExtractMxArrayField(matlabArray, fieldName, index);
    if(fieldPtr==NULL)
      throw ndlexceptions::MatlabInterfaceReadError(fieldName);
    return mxArrayToInt(fieldPtr);
  }
  double mxArrayExtractDoubleField(const mxArray* matlabArray, const string fieldName, const int index=0) const throw(ndlexceptions::MatlabInterfaceReadError&) 
  {
    mxArray* fieldPtr = mxArrayExtractMxArrayField(matlabArray, fieldName, index);
    if(fieldPtr==NULL)
      throw ndlexceptions::MatlabInterfaceReadError(fieldName);
    return mxArrayToDouble(fieldPtr);
  }
  string mxArrayExtractStringField(const mxArray* matlabArray, const string fieldName, const int index=0) const throw(ndlexceptions::MatlabInterfaceReadError&)
  {
    mxArray* fieldPtr = mxArrayExtractMxArrayField(matlabArray, fieldName, index);
    if(fieldPtr==NULL)
      throw ndlexceptions::MatlabInterfaceReadError(fieldName);
    return mxArrayToString(fieldPtr);
  }
  
  vector<int> mxArrayExtractVectorIntField(const mxArray* matlabArray, const string fieldName, const int index=0) const throw(ndlexceptions::MatlabInterfaceReadError&)
  {
    mxArray* fieldPtr = mxArrayExtractMxArrayField(matlabArray, fieldName, index);
    if(fieldPtr==NULL)
      throw ndlexceptions::MatlabInterfaceReadError(fieldName);
    return mxArrayToVectorInt(fieldPtr);
  }
  
  vector<unsigned int> mxArrayExtractVectorUintField(const mxArray* matlabArray, const string fieldName, int index=0) const throw(ndlexceptions::MatlabInterfaceReadError&)
  {
    mxArray* fieldPtr = mxArrayExtractMxArrayField(matlabArray, fieldName, index);
    if(fieldPtr==NULL)
      throw ndlexceptions::MatlabInterfaceReadError(fieldName);
    return mxArrayToVectorUint(fieldPtr);
  }
  bool mxArrayExtractBoolField(const mxArray* matlabArray, const string fieldName, const int index=0) const throw(ndlexceptions::MatlabInterfaceReadError&)
  {
    mxArray* fieldPtr = mxArrayExtractMxArrayField(matlabArray, fieldName, index);
    if(fieldPtr==NULL)
      throw ndlexceptions::MatlabInterfaceReadError(fieldName);
    return mxArrayToBool(fieldPtr);
  }
  
};

#else /* not _NDLMATLAB */

class CMatInterface {
 public:
 private:
};
#endif
#endif /* not CMATLAB_H*/
