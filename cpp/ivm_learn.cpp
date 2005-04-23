#include <fstream>
#include "CKern.h"
#include "CMatrix.h"
#include "CIvm.h"

int VERBOSITY=1;
int readSvmlDataFile(CMatrix& X, CMatrix& y, const string fileName);
void helpInfo();
void wait_any_key();
int main(int argc, char* argv[])
{
  string type="c";
  int maxiter=100;
  double tol=1e-6;
  double gamma=-1;
  int kernelType=0;
  double rbfInvWidth=1.0;
  string noiseType;
  int activeSetSize = 100;
  int i;
  for(i=1; (i<argc) && ((argv[i])[0] == '-');i++) {
    switch ((argv[i])[1]) 
      { 
      case '?': 
	helpInfo(); 
	exit(0);
      case 'd':
	i++;
	activeSetSize = atol(argv[i]);
	break;
      case 'z': 
	i++; 
	type = argv[i]; 
	break;
      case 'v': 
	i++; 
	VERBOSITY=atol(argv[i]); 
	break;
      case '#': 
	i++; 
	maxiter=atol(argv[i]); 
	break;
      case 'p': 
	i++; 
	gamma=atof(argv[i]); 
	break;
      case 't': 
	i++; 
	kernelType=atol(argv[i]); 
	break;
      case 'g': 
	i++; 
	rbfInvWidth=atof(argv[i]); 
	break;
      default: 
	cout << "Unrecognised option " << argv[i] << endl;
	helpInfo();
	exit(0);
      }
  }
  if(i>=argc) {
    cout << endl << "Not enough input parameters!" << endl << endl;
    wait_any_key();
    helpInfo();
    exit(0);
  }
  string trainData=argv[i];
  string modelfile;
  if((i+1)<argc) 
    modelfile=argv[i+1];
  if(type=="c")
    noiseType = "probit";
  else if(type=="r") 
    noiseType = "gaussian";
  else 
    {
      cout << "Unknown type " << type << ": Valid types are 'c' (classification), 'r' regession";
      wait_any_key();
      helpInfo();
      exit(0);
    }    
  if(gamma>1) 
    {
      cout << endl << "The fraction of unlabeled examples to classify as positives must" << endl;
      cout << "be less than 1.0 !!!" << endl << endl;
      wait_any_key();
      helpInfo();
      exit(0);
    }
  string dataFile="";
  CMatrix X;
  CMatrix y;
  readSvmlDataFile(X, y, trainData);

  CProbitNoise noise(y);

  
  CCmpndKern kern(X);
  switch(kernelType)
    {
    case 0:
      kern.addKern(new CLinKern(X));
      break;
    case 2:
      kern.addKern(new CRbfKern(X));
      // set rbf width to 2*gamma;
      break;
    }
  kern.addKern(new CBiasKern(X));
  kern.addKern(new CWhiteKern(X));
  
  int selectionCriterion = CIvm::ENTROPY; // entropy selection
  CIvm model(X, y, kern, noise, selectionCriterion, activeSetSize, VERBOSITY);
  model.selectPoints();
  //model.optimise(4, 100, 20);
  model.writeMatlabFile("testIvm.mat", "ivmInfo");
  model.kern.updateMatlabFile("testIvm.mat", "kern");
  model.noise.updateMatlabFile("testIvm.mat", "noise");
  X.updateMatlabFile("testIvm.mat", "X");
  y.updateMatlabFile("testIvm.mat", "y");

}
void helpInfo()
{
  cout << endl << "IVM Code: Version " <<endl;
  cout << "   usage: ivm_learn [options] example_file model_file" << endl << endl;
  cout << "Arguments:" << endl;
  cout << "         example_file-> file with training data" << endl;
  cout << "         model_file  -> file to store learned decision rule" << endl;
  cout << endl;
  cout << "General options:" << endl;
  cout << "         -?          -> this help" << endl;
  cout << "         -v [0..3]   -> verbosity level (default 1)" << endl;
  cout << "Learning options:" << endl;
  cout << "         -z {c,r}  -> select between classification and regression." << endl;
  cout << "Transduction options:" << endl;
  cout << "         -p [0..1] -> prior probability of an unlabelled data point" <<endl;
  cout << "                      being from the positive class (default is the ratio" << endl;
  cout << "                      of positive and negative examples in the training data)." << endl;
  cout << "Kernel options:" << endl;
  cout << "         -t int    -> type of kernel function:" << endl;
  cout << "                      0: linear (default)" << endl;
  cout << "                      2: radial basis function exp(-gamma ||a-b||^2)" << endl;
  cout << "         -g float  -> parameter gamma in rbf kernel." << endl;
}    
void wait_any_key()
{
  cout << endl << "(more)" << endl;
  (void)getc(stdin);
}

int oldTest()
{
  int numData = 1000;
  int numFeatures = 6;
  int numProcess = 1;
  CMatrix X(numData, numFeatures);
  X.randn();
  X.writeMatlabFile("test.mat", "X");

  // create weight vector.
  CMatrix w(numFeatures, 1);
  w.randn();
  w.updateMatlabFile("test.mat", "w");
  
  // create targets
  CMatrix y(numData, numProcess);
  // set random noise
  y.randn(0.01, 0.0);
  // add Xw to it
  y.gemv(X, w, 1.0, 1.0, "n");

  // get the sign of it.
  y.sign();
  y.updateMatlabFile("test.mat", "y");
  
  
}

int readSvmlDataFile(CMatrix& X, CMatrix& y, const string fileName)
{
  if(VERBOSITY>0)
    cout << "Loading training data." << endl;
  ifstream in(fileName.c_str());
  if(!in.is_open())
    throw("File not found");

  string line;
  string token;
  bool featureRead=false;
  int numData=0;
  int maxFeat=0;
  while(getline(in, line))
    {
      if(line[0]=='#')
	continue;
      else
	{
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
    }
  if(VERBOSITY>1)
    {
      cout << "Data number of features: " << maxFeat << endl;
      cout << "Number of data: " << numData << endl;
    }
  X.resize(numData, maxFeat);
  y.resize(numData, 1);
  in.close();
  ifstream inToo(fileName.c_str());
  int pointNo=0;
  while(getline(inToo, line))
    {
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
		      
		      string featStr=token.substr(0, ind);
		      string featValStr=token.substr(ind+1, token.size()-ind);
		      int featNum = atoi(featStr.c_str());
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
  if(VERBOSITY>0)
    cout << "Training data loaded." << endl;
}
