#include "CClctrl.h"

CClctrl::CClctrl(int arc, char** arv) : argc(arc), argv(arv)
{
  argNo=1;
  fileFormat = 0; // svmlight format
  verbosity = 2; // second highest level (highest is 3)
  time_t seconds;
  time(&seconds);
  setSeed((unsigned long) seconds);

}

bool CClctrl::isCurrentArg(string shortName, string commandName)
{
  return (getCurrentArgument()==shortName || getCurrentArgument()==commandName);
}
void CClctrl::confirmCurrentArg(string commandName)
{
  if(getCurrentArgument()!=commandName)
    exitError("Unrecognised argument " + getCurrentArgument() + " for `" + getMode() + "' command.");
}
void CClctrl::unrecognisedFlag()
{
  exitError("Unrecognised flag: " + getCurrentArgument() + " under " + getMode() +  " command.");
}
void CClctrl::helpUsage(const string description, const int width, const int padding)
{
  ndlstrutil::wrapOutputText(cout, "Usage: " + description, width, padding);
  cout << endl;
  cout << endl;
}
void CClctrl::helpDescriptor(const string description, const int width, const int padding)
{
  ndlstrutil::wrapOutputText(cout, description, width, padding);
  cout << endl;
  cout << endl;
}
void CClctrl::helpArgument(const string flags, const string explanation, const int width, const int padding)
{
  string padStr="";
  for(int i=0; i<padding; i++)
    padStr += " ";
  cout << padStr << flags << endl;
  ndlstrutil::wrapOutputText(cout, explanation, width, padding+3);
  cout << endl;
  cout << endl;
}
void CClctrl::waitForSpace()
{
  cout << "Press enter for more." << endl;
  char ch;
  cin.get(ch);    
}
int CClctrl::readSvmlDataFile(CMatrix& X, CMatrix& y, const string fileName)
{
  ifstream in(fileName.c_str());
  if(!in.is_open()) throw ndlexceptions::FileReadError(fileName);
  
  string line;
  string token;
  bool featureRead=false;
  int numData=0;
  int maxFeat=0;
  while(getline(in, line))
  {
    featureRead=false;
    
    if(line[line.size()-1]=='\r')
      line.erase(line.size()-1);
    if(line[0]=='#')
      continue;
    numData++;
    int pos=0;
    while(pos<line.size())
    {
      token.erase();
      while(pos<line.size() && line[pos]!=' ')
      {
	token+=line[pos];
	pos++;
      }
      pos++;
      if(token.size()>0)
      {
	// deal with token.
	if(featureRead)
	{
	  int ind = token.find(':');
	  if(ind==std::string::npos || ind < 0)
	  {
	    ndlexceptions::StreamFormatError err("");
	    throw ndlexceptions::FileFormatError(fileName, err);
	  }
	  string featStr=token.substr(0, ind);
	  int featNum = atoi(featStr.c_str());
	  if(featNum>maxFeat)
	    maxFeat=featNum;
	}
	else
	{
	  featureRead=true;
	}
      }
    }
    
  }
  if(verbosity>1)
  {
    cout << "Data number of features: " << maxFeat << endl;
    cout << "Number of data: " << numData << endl;
  }
  X.resize(numData, maxFeat);
  X.zeros();
  y.resize(numData, 1);
  y.zeros();
  in.close();
  ifstream inToo(fileName.c_str());
  int pointNo=0;
  while(getline(inToo, line))
  {
    if(line[line.size()-1]=='\r')
      line.erase(line.size()-1);
    featureRead=false;
    if(line[0]=='#')
      continue;
    else
    {
      int pos=0;
      while(pos<line.size())
      {
	token.erase();
	while(pos<line.size() && line[pos]!=' ')
	{
	  token+=line[pos];
	  pos++;
	}
	pos++;
	if(token.size()>0)
	{
	  // deal with token.
	  if(featureRead)
	  {
	    int ind = token.find(':');		      
	    // TODO Check that : is in the string.
	    string featStr=token.substr(0, ind);
	    string featValStr=token.substr(ind+1, token.size()-ind);
	    int featNum = atoi(featStr.c_str());
	    if(featNum<1 || featNum>maxFeat || pointNo<0 || pointNo>=numData)
	    {
	      ndlexceptions::StreamFormatError err("");
	      throw ndlexceptions::FileFormatError(fileName, err);
	    }
		      
	    double featVal = atof(featValStr.c_str());
	    X.setVal(featVal, pointNo, featNum-1);
	  }
	  else
	  {
	    y.setVal(atof(token.c_str()), pointNo);
	    featureRead=true;
	  }
	}
      }
      
    }
    pointNo++;
  }
  // WVB ADDED
  return -1;
}

void CClctrl::readData(CMatrix& X, CMatrix& y, const string fileName)
{
  string m = getMode();
  setMode("file");
  if(verbosity>1)
    cout << "Loading data." << endl;
  switch(fileFormat)
  {
  case 0: /// svmlight file format.
    readSvmlDataFile(X, y, fileName);
    break;
  case 1: /// Matlab file format.
#ifdef _NDLMATLAB
    X.readMatlabFile(fileName, "X");
    y.readMatlabFile(fileName, "y");
#else
    throw ndlexceptions::MatlabInterfaceError("MATLAB not incorporated at compile time");
#endif
    break;
  default:
    exitError("Unrecognised file format number.");
    
  }
  if(verbosity>1)
    cout << "Data set loaded." << endl;
  setMode(m);
}
void CClctrl::exitNormal()
{
#ifdef _DEBUG
#ifdef _MSC_VER
  // For debugging under visual studio, to prevent window disappearing.
  waitForSpace();
#endif
#endif
  exit(0);
}
void CClctrl::exitError(const string error)
{
  cerr << error << endl << endl;
  waitForSpace();
  helpInfo();
  exit(1);
}

