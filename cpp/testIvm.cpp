#include "CKern.h"
#include "CMatrix.h"
#include "CIvm.h"

int main()
{
  int numData = 20;
  int numFeatures = 3;
  int numProcess = 1;
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
  CMatrix X(numFeatures, numData, x);
  X.trans();
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
  CIvm model(X, y, kern, noise, 0, 10);
  for(int iters=0; iters<10; iters++)
    {
      model.init();
      model.selectPoints();
      for(int i=0; i<model.activeSet.size(); i++)
	cout << model.activeSet[i] << endl;
      model.checkGradients();
      model.gdOptimise(0.01, 0.9, 0, 1000);
    }
  model.init();
  model.selectPoints();
  kern.display(cout);
}


