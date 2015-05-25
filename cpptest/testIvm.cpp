#include "CKern.h"
#include "CMatrix.h"
#include "CIvm.h"

int testGaussian();
int testProbit();
int testNcnm();

int main()
{
  int fail = 0;
  fail += testGaussian();
  fail += testProbit();
  fail += testNcnm();
}

int testGaussian()
{
  int fail = 0;
  CMatrix X;
  X.readMatlabFile("testGaussian.mat", "X");
  CMatrix y;
  y.readMatlabFile("testGaussian.mat", "y");

  CCmpndKern kernInit(X);
  kernInit.addKern(new CRbfKern(X));
  kernInit.addKern(new CLinKern(X));
  kernInit.addKern(new CBiasKern(X));
  kernInit.addKern(new CWhiteKern(X));
  CGaussianNoise noiseInit(y);
  
  CIvm modelInit(X, y, kernInit, noiseInit, CIvm::ENTROPY, 50);
  modelInit.selectPoints();
   
  CCmpndKern kern(X);
  kern.readMatlabFile("testGaussian.mat", "kernInit");
  CGaussianNoise noise(y);
  noise.readMatlabFile("testGaussian.mat", "noiseInit");

  CIvm model(X, y, kern, noise, "testGaussian.mat", "ivmInfoInit", 0);
  if(model.equals(modelInit))
    cout << model.getNoiseName() << " Noise IVM passed." << endl;
  else
    {
      cout << "FAILURE: " << model.getNoiseName() << " Noise IVM." << endl;
      fail++;
    }
  return fail;
}
int testProbit()
{
  int fail = 0;
  CMatrix X;
  X.readMatlabFile("testProbit.mat", "X");
  CMatrix y;
  y.readMatlabFile("testProbit.mat", "y");

  CCmpndKern kernInit(X);
  kernInit.addKern(new CRbfKern(X));
  kernInit.addKern(new CLinKern(X));
  kernInit.addKern(new CBiasKern(X));
  kernInit.addKern(new CWhiteKern(X));
  CProbitNoise noiseInit(y);
  
  CIvm modelInit(X, y, kernInit, noiseInit, CIvm::ENTROPY, 50);
  modelInit.selectPoints();
   
  CCmpndKern kern(X);
  kern.readMatlabFile("testProbit.mat", "kernInit");
  CProbitNoise noise(y);
  noise.readMatlabFile("testProbit.mat", "noiseInit");

  CIvm model(X, y, kern, noise, "testProbit.mat", "ivmInfoInit", 0);
  if(model.equals(modelInit))
    cout << model.getNoiseName() << " Noise IVM passed." << endl;
  else
    {
      cout << "FAILURE: " << model.getNoiseName() << " Noise IVM." << endl;
      fail++;
    }
  return fail;
}

int testNcnm()
{
  int fail = 0;
  CMatrix X;
  X.readMatlabFile("testNcnm.mat", "X");
  CMatrix y;
  y.readMatlabFile("testNcnm.mat", "y");

  CCmpndKern kernInit(X);
  kernInit.addKern(new CRbfKern(X));
  kernInit.addKern(new CLinKern(X));
  kernInit.addKern(new CBiasKern(X));
  kernInit.addKern(new CWhiteKern(X));

  // Add L1 prior to the kernel.
  CGammaDist* prior = new CGammaDist();
  prior->setParam(1.0, 0);
  prior->setParam(1.0, 1);
  kernInit.addPrior(prior, 0);
  kernInit.addPrior(prior, 2);
  kernInit.addPrior(prior, 3);
  kernInit.addPrior(prior, 4);
  CNcnmNoise noiseInit(y);
  
  CIvm modelInit(X, y, kernInit, noiseInit, CIvm::ENTROPY, 50);
  modelInit.selectPoints();
  modelInit.checkGradients(); 
  CCmpndKern kern(X);
  kern.readMatlabFile("testNcnm.mat", "kernInit");
  CNcnmNoise noise(y);
  noise.readMatlabFile("testNcnm.mat", "noiseInit");
  
  CIvm model(X, y, kern, noise, "testNcnm.mat", "ivmInfoInit", 0);
  if(model.equals(modelInit))
    cout << model.getNoiseName() << " Noise IVM passed." << endl;
  else
    {
      cout << "FAILURE: " << model.getNoiseName() << " Noise IVM." << endl;
      fail++;
    }
  return fail;
}


