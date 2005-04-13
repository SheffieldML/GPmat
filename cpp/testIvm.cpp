#include "CKern.h"
#include "CMatrix.h"
#include "CIvm.h"

int main()
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
  y.updateMatlabFile("test.mat", "y");
  
  CGaussianNoise noise(y);


  CCmpndKern kern(X);
  kern.addKern(new CRbfKern(X));
  kern.addKern(new CLinKern(X));
  kern.addKern(new CBiasKern(X));
  kern.addKern(new CWhiteKern(X));
  cout << "Kernel:" << endl;
  kern.display(cout);
  int verbosity = 2;
  int activeSetSize = 100;
  int selectionCriterion = 0; // entropy selection
  CIvm model(X, y, kern, noise, selectionCriterion, activeSetSize, verbosity);
  model.optimise(4, 100, 20);
}


