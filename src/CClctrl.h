
#ifndef CCLCTRL_H
#define CCLCTRL_H
#include <iostream>
#include <ctime>
#include <cstring>
#include "ndlstrutil.h"
#include "ndlutil.h"
#include "ndlfortran.h"
#include "ndlexceptions.h"
#include "CMatrix.h"

// Command line control class header.
using namespace std;
class CClctrl {
 public:
  // Constructor given the input arguments.
  CClctrl(int argc, char** argv);
  virtual ~CClctrl() {}
  // Check if the current argument is one of the given flags
  bool isCurrentArg(string shortName, string longName);
  // Check that the current argument is valid.
  void confirmCurrentArg(string argument);
  void unrecognisedFlag();
  // helper function for formatting help information.
  void helpUsage(const string description, const int width=79, const int padding=0);
  // helper function for formating help information.
  void helpDescriptor(const string description, const int width=79, const int padding=5);
  // helper function for formating help information.
  void helpArgument(const string flags, const string explanation, const int width=79, const int padding=5);
  // pauses the screen until user presses space.
  void waitForSpace();
  // Read in a data file in Thorsten Joachim's SVM light format.
  int readSvmlDataFile(CMatrix& X, CMatrix& y, const string fileName);
  // Read in a data file.
  void readData(CMatrix& X, CMatrix& y, const string fileName);
  // exit normally.
  void exitNormal();
  // Exit with an error.
  void exitError(const string error);

  // Virtual functions for help commands.
  virtual void helpInfo()=0;
  virtual void helpHeader()=0;
  
  // true if there are still unread flags from the command line.
  bool isFlags() const
    {
      return flags && getCurrentArgumentNo()<argc;
    }
  void setFlags(bool val) 
    {
      flags = val;
    }

  // get and set the verbosity level.
  int getVerbosity() const
    {
      return verbosity;
    }
  void setVerbosity(int val)
    {
      verbosity = val;
    }

  // get and set the seed for the random number generator.
  unsigned long getSeed() const
    {
      return seed;
    }
  void setSeed(unsigned long val)
    {
      ndlutil::init_genrand(val);
      seed = val;
    }

  // get and set the file format for input files.
  int getFileFormat() const
    {
      return fileFormat;
    }
  void setFileFormat(int val)
    {
      fileFormat = val;
    }

  // Manipulate the current command line argument number
  void incrementArgument()
    {
      argNo++;
    }
  int getCurrentArgumentNo() const
    {
      return argNo;
    }
  void setCurrentArgumentNo(int val)
    {
      argNo=val;
    }

  // Recover information from the command line arguments.
  string getCurrentArgument() const
    {
      BOUNDCHECK(argNo<argc && argNo>=0);
      return argv[argNo];
    }
  int getIntFromCurrentArgument() const
    {
      BOUNDCHECK(argNo<argc && argNo>=0);
      return atol(argv[argNo]);
    }
  double getDoubleFromCurrentArgument() const
    {
      BOUNDCHECK(argNo<argc && argNo>=0);
      return atof(argv[argNo]);
    }
  bool getBoolFromCurrentArgument() const
    {
      BOUNDCHECK(argNo<argc && argNo>=0);
      // Assume it is a string first.
      string arg = getStringFromCurrentArgument();
      if(arg=="true" || arg=="1")
	return true;
      else if(arg=="false" || arg=="0")
	return false;
      else 
	throw ndlexceptions::CommandLineError("Current argument is not boolean.");
    }
  string getStringFromCurrentArgument() const
    {
      BOUNDCHECK(argNo<argc && argNo>=0);
      return argv[argNo];
    }
  double getCurrentArgumentLength() const
    {
      return strlen(argv[argNo]);
    }
  // test if the argument is a flag (i.e. starts with -)
  bool isCurrentArgumentFlag() const
    {
      if(argv[argNo][0]=='-')
	return true;
      else
	return false;
    }
  // manipulate the mode --- typically for determining what the help file output will be 
  void setMode(string val) 
    {
      mode = val;
    }
  string getMode() const
    {
      return mode;
    }

 private:
  bool flags;
  unsigned long seed;
  int verbosity;
  int fileFormat;
  int argNo;
  string mode;
 protected:
  int argc; 
  char** argv; 
  
};

#endif

