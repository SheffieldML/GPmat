#include "CNoise.h"


int main()
{
  int numProcess = 3;
  int numData = 1000;
  CMatrix target(numData, numProcess);
  target.randn();
  CGaussianNoise noise(target);
  CMatrix params(1, noise.getNParams());
  noise.getTransParams(params);
  cout << "Original parameters: " << endl << params << endl;
  params.randn();
  double varVal = randn();
  // make sure the variance is positive
  params.setVals(varVal*varVal, 0, params.getCols()-1);
  noise.setTransParams(params);
  cout << "New parameters: " << endl << params << endl;

    
  CMatrix mu(numData, numProcess);
  mu.randn();
  CMatrix varsigma(numData, numProcess);
  varsigma.randn();
  varsigma*=varsigma;
  //CMatrix y(numData, numProcess);
  CMatrix y = noise.out(mu, varsigma);
  y.randn();
  cout << "mu values" << endl << mu << endl;
  cout << "varsigma values" << endl << varsigma << endl;
  cout << "y values " << endl << y << endl;

  double epsilon = 1e-6;
  double Lplus;
  double Lminus;
  CMatrix diffGradParam(1, params.getCols());
  CMatrix origParams(1, params.getCols());
  noise.getTransParams(origParams);
  for(int j=0; j<params.getCols(); j++)
    {
      params.setVals(origParams.getVals(j) + epsilon, j);
      noise.setTransParams(params);
      Lplus = noise.logLikelihood(mu, varsigma, y);
      params.setVals(origParams.getVals(j) - epsilon, j);
      noise.setTransParams(params);
      Lminus = noise.logLikelihood(mu, varsigma, y);
      diffGradParam.setVals(0.5*(Lplus - Lminus)/epsilon, j);
      params.setVals(origParams.getVals(j), j);
    }
  
  cout << "Numerical differences:" << endl << diffGradParam << endl;
  CMatrix analyticalGradParam(1, params.getCols());
  noise.getGradTransParams(analyticalGradParam, mu, varsigma, y);
  cout << "Analytic gradients:" << endl << analyticalGradParam << endl;
  CMatrix diffParam(params.getRows(), params.getCols());
  diffParam.deepCopy(analyticalGradParam);
  diffParam-=diffGradParam;
  cout << "Maximum param difference: " << diffParam.max() << endl;

  // Check mu gradients.
  CMatrix diffGradMu(mu.getRows(), mu.getCols());
  CMatrix origMu(mu.getRows(), mu.getCols());
  origMu.deepCopy(mu);
  for(int j=0; j<mu.getNumElements(); j++)
    {
      mu.setVals(origMu.getVals(j) + epsilon, j);
      Lplus = noise.logLikelihood(mu, varsigma, y);
      mu.setVals(origMu.getVals(j) - epsilon, j);
      Lminus = noise.logLikelihood(mu, varsigma, y);
      diffGradMu.setVals(0.5*(Lplus - Lminus)/epsilon, j);
      mu.setVals(origMu.getVals(j), j);
    }
  // Check varsigma gradients.
  CMatrix diffGradVarsigma(varsigma.getRows(), varsigma.getCols());
  CMatrix origVarsigma(varsigma.getRows(), varsigma.getCols());
  origVarsigma.deepCopy(varsigma);
  for(int j=0; j<varsigma.getNumElements(); j++)
    {
      varsigma.setVals(origVarsigma.getVals(j) + epsilon, j);
      Lplus = noise.logLikelihood(mu, varsigma, y);
      varsigma.setVals(origVarsigma.getVals(j) - epsilon, j);
      Lminus = noise.logLikelihood(mu, varsigma, y);
      diffGradVarsigma.setVals(0.5*(Lplus - Lminus)/epsilon, j);
      varsigma.setVals(origVarsigma.getVals(j), j);
    }
  
  CMatrix analyticalGradMu(mu.getRows(), mu.getCols());
  CMatrix analyticalGradVarsigma(varsigma.getRows(), varsigma.getCols());
  noise.getGradInputs(analyticalGradMu, analyticalGradVarsigma, mu, varsigma, y);
  cout << endl;
  cout << "Mu numerical differences " << endl << diffGradMu << endl;
  cout << "Mu analytical gradient " << endl << analyticalGradMu << endl;

  CMatrix diffMu(mu.getRows(), mu.getCols());
  diffMu.deepCopy(analyticalGradMu);
  diffMu-=diffGradMu;
  cout << "Maximum mu difference: " << diffMu.max() << endl;
  cout << endl;
  cout << "Varsigma numerical differences " << endl << diffGradVarsigma << endl;
  cout << "Varsigma analytical gradient " << endl << analyticalGradVarsigma << endl;
  CMatrix diffVarsigma(varsigma.getRows(), varsigma.getCols());
  diffVarsigma.deepCopy(analyticalGradVarsigma);
  diffVarsigma-=diffGradVarsigma;
  cout << "Maximum varsigma difference: " << diffVarsigma.max() << endl;


}
