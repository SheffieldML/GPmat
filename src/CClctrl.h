#ifndef CCLCTRL_H
#define CCLCTRL_H

#include "ndlstrutil.h"
#include "ndlexceptions.h"
#include "CMatrix.h"
// Command line control class header.
using namespace std;
class CClctrl {
 public:
  CClctrl(){}
  CClctrl(int argc, char** argv);
  void confirmCurrentArg(string argument);
  void unrecognisedFlag();
  void helpUsage(const string description, const int width=79, const int padding=0);
  void helpDescriptor(const string description, const int width=79, const int padding=5);
  void helpArgument(const string flags, const string explanation, const int width=79, const int padding=5);
  void waitForSpace();
  int readSvmlDataFile(CMatrix& X, CMatrix& y, const string fileName);
  void readData(CMatrix& X, CMatrix& y, const string fileName);
  void exitError(const string error);

  virtual void helpInfo()=0;
  virtual void helpHeader()=0;

  bool isFlags() const
    {
      return flags && getCurrentArgumentNo()<argc;
    }
  void setFlags(bool val) 
    {
      flags = val;
    }
  int getVerbosity() const
    {
      return verbosity;
    }
  void setVerbosity(int val)
    {
      verbosity = val;
    }
  int getFileFormat() const
    {
      return fileFormat;
    }
  void setFileFormat(int val)
    {
      fileFormat = val;
    }
  void incrementArgument()
    {
      argNo++;
    }
  string getCurrentArgument() const
    {
      assert(argNo<argc && argNo>=0);
      return argv[argNo];
    }
  int getIntFromCurrentArgument() const
    {
      assert(argNo<argc && argNo>=0);
      return atol(argv[argNo]);
    }
  double getDoubleFromCurrentArgument() const
    {
      assert(argNo<argc && argNo>=0);
      return atof(argv[argNo]);
    }
  bool getBoolFromCurrentArgument() const
    {
      assert(argNo<argc && argNo>=0);
      int val = atol(argv[argNo]);
      if(val==0)
	return false;
      else if(val==1)
	return true;
      else 
	throw ndlexceptions::CommandLineError("Current argument is not boolean");
    }
	
  string getStringFromCurrentArgument() const
    {
      assert(argNo<argc && argNo>=0);
      return argv[argNo];
    }
  double getCurrentArgumentLength() const
    {
      return strlen(argv[argNo]);
    }
  int getCurrentArgumentNo() const
    {
      return argNo;
    }
  void setCurrentArgumentNo(int val)
    {
      argNo=val;
    }
  bool isCurrentArgumentFlag() const
    {
      if(argv[argNo][0]=='-')
	return true;
      else
	return false;
    }
  string getFlagText() const
    {
      return argv[argNo]+1;
    }
  void setMode(string val) 
    {
      mode = val;
    }
  string getMode() const
    {
      return mode;
    }
  void setArgv(char** arg)
    {
      argv = arg;
    }
  void setArgc(int arc)
    {
      argc = arc;
    }
 private:
  bool flags;
  int verbosity;
  int fileFormat;
  int argNo;
  string mode;
 protected:
  int argc; 
  char** argv; 
  
};

#endif
