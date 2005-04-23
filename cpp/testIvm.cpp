#include "CKern.h"
#include "CMatrix.h"
#include "CIvm.h"

int testGaussian();
int testProbit();

int main()
{
  int fail = 0;
  fail += testGaussian();
  fail += testProbit();
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
    cout << "Gaussian Noise passed." << endl;
  else
    {
      cout << "FAILURE: Gaussian noise." << endl;
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
    cout << "Probit Noise passed." << endl;
  else
    {
      cout << "FAILURE: Probit noise." << endl;
      fail++;
    }
  return fail;
}


