#include "CKern.h"
using namespace std;

int testType(const string kernelType);
int testKern(CKern* kern, CKern* kern2, const string fileName);

int main()
{
  int fail=0;
  char* testKern = "bias";
  fail += testType(testKern);
   fail += testType("white");
   fail += testType("lin");
   fail += testType("rbf");
   fail += testType("cmpnd");
  cout << "Number of failures: " << fail << "." << endl;
}

int testType(const string kernelType)
{
  string fileName = kernelType + "Test.mat";

  CMatrix X;
  X.readMatlabFile(fileName, "X");

  CKern* kern;
  CKern* kern2;
 
  if(kernelType=="lin")
    {
      kern = new CLinKern(X);
      kern2 = new CLinKern(X);
    }
  else if(kernelType=="rbf")
    {
      kern = new CRbfKern(X);
      kern2 = new CRbfKern(X);
    }
  else if(kernelType=="white")
    {
      kern = new CWhiteKern(X);
      kern2 = new CWhiteKern(X);
    }
  else if(kernelType=="bias")
    {
      kern = new CBiasKern(X);
      kern2 = new CBiasKern(X);
    }
  else if(kernelType=="cmpnd")
    {
      kern = new CCmpndKern(X);
      kern->addKern(new CRbfKern(X));
      kern->addKern(new CLinKern(X));
      kern->addKern(new CBiasKern(X));
      kern->addKern(new CWhiteKern(X));
      kern2 = new CCmpndKern(X);

    }
  int fail = testKern(kern, kern2, fileName);
  delete kern;
  delete kern2;
  return fail;
}
int testKern(CKern* kern, CKern* kern2, const string fileName)
{
  
  int fail = 0;
  CMatrix params;
  params.readMatlabFile(fileName, "params");
  CMatrix X;
  X.readMatlabFile(fileName, "X");
  CMatrix X2;
  X2.readMatlabFile(fileName, "X2");
  kern->setTransParams(params);
  kern2->readMatlabFile(fileName, "kern2");
  if(kern2->equals(*kern))
    cout << kern->getName() << " Initial Kernel matches." << endl;
  else
    {
      cout << "FAILURE: " << kern->getName() << " Initial Kernel." << endl;
      fail++;
    }
  CMatrix K1(X.getRows(), X.getRows());
  kern->compute(K1, X);
  CMatrix K2;
  K2.readMatlabFile(fileName, "K2");
  if(K1.equals(K2))
    cout << kern->getName() << " full compute matches." << endl;
  else
    { 
      cout << "FAILURE: " << kern->getName() << " full compute." << endl;
      fail++;
    }
   CMatrix K3(X.getRows(), X2.getRows());
   kern->compute(K3, X, X2);
   CMatrix K4;
   K4.readMatlabFile(fileName, "K4");
   if(K3.equals(K4))
     cout << kern->getName() << " double compute matches." << endl;
   else
     { 
       cout << "FAILURE: " << kern->getName() << " double compute." << endl;
       fail++;
     }
   CMatrix k1(X.getRows(), 1);
   kern->diagCompute(k1, X);
   CMatrix k2;
   k2.readMatlabFile(fileName, "k2");
   if(k1.equals(k2))
     cout << kern->getName() << " diag compute matches." << endl;
   else
     {
       cout << "FAILURE: " << kern->getName() << " diag compute." << endl;
       fail++;
     }
   CMatrix covGrad;
   covGrad.readMatlabFile(fileName, "covGrad");
   covGrad.setSymmetric(true);
   CMatrix g1(1, kern->getNumParams());
   kern->getGradTransParams(g1, X, covGrad);
   CMatrix g2;
   g2.readMatlabFile(fileName, "g2");
   if(g1.equals(g2))
     cout << kern->getName() << " parameter gradient matches." << endl;
   else
     {
       cout << "FAILURE: " << kern->getName() << " parameter gradient." << endl;
       fail++;
     }
   
   kern->writeMatlabFile("crap.mat", "writtenKern");
   kern2->readMatlabFile("crap.mat", "writtenKern");
   if(kern->equals(*kern2))
     cout << "Written kernel matches read in kernel. Read and write to matlab passes." << endl;
   else
     {
       cout << "FAILURE: Read in kernel does not match written out kernel." << endl;
       fail++;
     }
   return fail;
}
  

/*
  int nrows = 3;
  int ncols = 20;
  double x[60] = {-0.4326,    0.2944,   -1.6041,
		  -1.6656,   -1.3362,   0.2573,
		  0.1253,    0.7143,   -1.0565,
		  0.2877,    1.6236,    1.4151,
		  -1.1465,   -0.6918,   -0.8051,
		  1.1909,    0.8580,    0.5287,
		  1.1892,    1.2540,    0.2193,
		  -0.0376,   -1.5937,   -0.9219,
		  0.3273,   -1.4410,   -2.1707,
		  0.1746,    0.5711,   -0.0592,
		  -0.1867,   -0.3999,   -1.0106,
		  0.7258,    0.6900,    0.6145,
		  -0.5883,    0.8156,    0.5077,
		  2.1832,    0.7119,    1.6924,
		  -0.1364,    1.2902,    0.5913,
		  0.1139,    0.6686,   -0.6436,
		  1.0668,    1.1908,    0.3803,
		  0.0593,   -1.2025,   -1.0091,
		  -0.0956,   -0.0198,   -0.0195,
		  -0.8323,   -0.1567,   -0.0482};
  CMatrix X(nrows, ncols, x);
  X.trans();
  X.updateMatlabFile("testKern.mat", "X");
  CWhiteKern kernWhite(X);
  CMatrix K = kernWhite.compute(X);
  //  cout << "Kernel: " << endl << K << endl;
  CBiasKern kernBias(X);
  CMatrix K2 = kernBias.compute(X);
  K += K2;
  //  cout << "Kernel: " << endl << K << endl;
  CLinKern kernLin(X);
  CMatrix K3 = kernLin.compute(X);
  K += K3;
  // cout << "Kernel: " << endl << K << endl;
  kernLin.display(cout);// << kernLin << endl;//  << kernWhite << kernBias << endl;
  kernWhite.display(cout);
  kernBias.display(cout);
   CMatrix Kinv(K.getRows(), K.getCols());
  Kinv.deepCopy(K);
  Kinv.pdinv();
  cout << Kinv;
  CCmpndKern kernCmpnd(X);
  kernCmpnd.addKern(new CLinKern(X));
  kernCmpnd.addKern(new CWhiteKern(X));
  kernCmpnd.addKern(new CBiasKern(X));
  kernCmpnd.display(cout);//cout << "Kernel: "<< endl << kernCmpnd;
  kernBias.updateMatlabFile("testKern.mat", "kernBias");
  kernCmpnd.updateMatlabFile("testKern.mat", "kernCmpnd");
  CMatrix K4 = kernCmpnd.compute(X);
  K4.updateMatlabFile("testKern.mat", "K");

  K-=K4;
  cout << endl << K;
  K4.pdinv();
  cout << endl << K4;

  cout << "First row and first column of knerCmpnd is  ";
  cout << kernCmpnd.computeElement(X, 0, X, 0) << endl;
  CCmpndKern testCmpnd(X);
  testCmpnd.readMatlabFile("testKern.mat", "kernCmpnd");
  testCmpnd.writeMatlabFile("testKern2.mat", "kernCmpnd");
*/
