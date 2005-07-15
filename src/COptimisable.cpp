#include "COptimisable.h"

void COptimisable::checkGradients()
{
  double change = ndlutil::GRADCHANGE;
  double origParam = 0.0;
  double objectivePlus = 0.0;
  double objectiveMinus = 0.0;
  
  CMatrix analyticGrad(1, getOptNumParams());
  CMatrix numericalDiff(1, getOptNumParams());
  CMatrix diffNumericalAnalytic(1, getOptNumParams());
  CMatrix params(1, getOptNumParams());
  CMatrix origParams(1, getOptNumParams());
  
  getOptParams(params);
  origParams.deepCopy(params);
  computeObjectiveGradParams(analyticGrad);
  for(int j=0; j<getOptNumParams(); j++)
    {
      origParam = origParams.getVal(j);
      params.setVal(origParam + change, j);
      setOptParams(params);
      objectivePlus = computeObjectiveVal();
      params.setVal(origParam - change, j);
      setOptParams(params);
      objectiveMinus = computeObjectiveVal();
      numericalDiff.setVal(0.5*(objectivePlus - objectiveMinus)/change, j);
      params.setVal(origParams.getVal(j), j);
    }
  
  cout << "Numerical differences:" << endl << numericalDiff << endl;
  cout << "Analytic gradients:" << endl << analyticGrad << endl;
  diffNumericalAnalytic.deepCopy(analyticGrad);
  diffNumericalAnalytic-=numericalDiff;
  cout << "Differences: " << endl << diffNumericalAnalytic << endl;
}

void COptimisable::gdOptimise(double learnRate, double momentum, int display, int maxIters, const double objectiveTol, const double paramTol)
{
  int nParams = getOptNumParams();
  double objectiveVal = 0.0;
  double oldObjective = 0.0;
  double diffObjective = 0.0;
  double diffParam = 0.0;
  CMatrix params(1, getOptNumParams());
  CMatrix oldParams(1, getOptNumParams());
  CMatrix gradParams(1, getOptNumParams());
  CMatrix changeParams(1, getOptNumParams());
  changeParams.zeros();
  getOptParams(params);
  if(evalFunc)
    objectiveVal = computeObjectiveVal();
  for(int iter=0; iter<maxIters; iter++)     
    {
      oldParams.deepCopy(params);
      computeObjectiveGradParams(gradParams);
      if(momentum>0)
	{
	  changeParams.axpy(gradParams, -learnRate/momentum);
	  params.axpy(changeParams, momentum);
	  changeParams.scale(momentum);
	}
      else
	{
	  params.axpy(gradParams, -learnRate);	    
	}
      setOptParams(params);
      if(evalFunc)
	{
	  oldObjective = objectiveVal;
	  objectiveVal = computeObjectiveVal();
	  diffObjective = abs(objectiveVal-oldObjective);
	}
      if(getVerbosity()>2)
	{
	  cout << "Iteration: " << iter << ", objective function: " << objectiveVal << endl;
	}
      diffParam = params.maxAbsDiff(oldParams);
      if(diffObjective<objectiveTol && diffParam<parameterTol)
	{
	  cout << "Param difference: " << diffParam << endl;
	  cout << "Objective difference: " << diffObjective << endl;
	  cout << "Converged .." << endl;
	  break;
	}
      
    } 
  cout << "Parameters: " << endl;
  cout << params << endl;
  /*       if(iter>=maxIters) */
  /* 	cout << "Maximum iterations exceed in gdOptimise(), objective change, " << diffObjective << ", parameter change, " << diffParam << endl; */
}
void COptimisable::gdPullbackOptimise(double learnRate, int display, int maxIters, const double objectiveTol, const double paramTol)
{
  int nParams = getOptNumParams();
  double objectiveVal = 0.0;
  double oldObjective = 0.0;
  double diffObjective = 0.0;
  double diffParam = 0.0;
  CMatrix params(1, getOptNumParams());
  CMatrix oldParams(1, getOptNumParams());
  CMatrix gradParams(1, getOptNumParams());
  CMatrix changeParams(1, getOptNumParams());
  changeParams.zeros();
  getOptParams(params);
  objectiveVal = computeObjectiveVal();
  for(int iter=0; iter<maxIters; iter++)     
    {
      while(1)
	{
	  oldParams.deepCopy(params);
	  computeObjectiveGradParams(gradParams);
	  params.axpy(gradParams, -learnRate);	    
	  setOptParams(params);
	  oldObjective = objectiveVal;
	  objectiveVal = computeObjectiveVal();
	  diffObjective = oldObjective-objectiveVal;
	  if(diffObjective<0)
	    {
	      params.deepCopy(oldParams);
	      learnRate = learnRate/2;
	      objectiveVal=oldObjective;
	      setOptParams(params);
	      if(getVerbosity()>2)
		cout << "Back tracking, learning rate: " << learnRate << endl;
	    }
	  else
	    {
	      learnRate = learnRate*1.1;
	      break;
	    }
	}
      if(getVerbosity()>2)
	{
	  cout << "Iteration: " << iter << ", objective function: " << objectiveVal << endl;
	}
      diffParam = params.maxAbsDiff(oldParams);
      if(diffObjective<objectiveTol && diffParam<parameterTol)
	{
	  cout << "Param difference: " << diffParam << endl;
	  cout << "Objective difference: " << diffObjective << endl;
	  cout << "Converged .." << endl;
	  break;
	}
      
    } 
  cout << "Parameters: " << endl;
  cout << params << endl;
  /*       if(iter>=maxIters) */
  /* 	cout << "Maximum iterations exceed in gdOptimise(), objective change, " << diffObjective << ", parameter change, " << diffParam << endl; */
}
double COptimisable::oneDObjectiveVal(const double val)
{
  assert(paramStoreOne.getRows()==1);
  assert(paramStoreTwo.getRows()==1);
  assert(paramStoreOne.getCols()==getOptNumParams());
  assert(paramStoreTwo.getCols()==getOptNumParams());
  getOptParams(paramStoreOne);
  paramStoreTwo.deepCopy(paramStoreOne);
  // add val times direction to the parameters.
  paramStoreTwo.axpy(direction, val);
  setOptParams(paramStoreTwo);
  double objective = computeObjectiveVal();
  setOptParams(paramStoreOne);
  return objective;
}
void COptimisable::scgOptimise(int maxIters, const double objectiveTol, const double paramTol)
{

  // taken from the paper by Martin Moller: "A scaled conjugate gradient algorithm for fast supervised learning".
  int nParams = getOptNumParams();
  double objectiveVal = 0.0;
  double oldObjective = 0.0;
  double diffObjective = 0.0;
  double diffParam = 0.0;
  
  bool success = true;
  double beta;
  double Delta;
  double lambda;
  double lambdaBar;
  double sigma;
  double rr;

  CMatrix w(1, getOptNumParams());
  CMatrix wPlus(1, getOptNumParams());

  CMatrix r(1, getOptNumParams());
  CMatrix p(1, getOptNumParams());
  CMatrix rp(1, getOptNumParams());
  CMatrix s(1, getOptNumParams());
  getOptParams(w);

  const double m_step = 1.0e-4;
  const double m_reg = 1.0;
  double oldObj = computeObjectiveVal();
  double obj = oldObj;
  double mu = 0.0;
  double theta = 0.0;
  double sigmaInv = 0.0;
  double delta = 0.0;
  double alpha = 0.0;
  double newObj = 0.0;
  double gamma = 0.0;
  int j = 1;					// j counts number of iterations.
  
  // 1
  lambda = m_reg; // lambda is the scale??
  lambdaBar = 0.0;
  computeObjectiveGradParams(r);
  r.negate();
  p.deepCopy(r);

  for(int k=1; k<=maxIters; k++)
    {
      // should check that p is not infinite above
      double normp = p.normRow(0);
      double normp2 = normp*normp;
      // 2
      if(success)
	{
	  // can get a divide by zero here if pp is too small.
	  
	  sigma = m_step/normp;
	  wPlus.deepCopy(w);
	  wPlus.axpy(p, sigma);
	  setOptParams(wPlus);
	  computeObjectiveGradParams(s);
	  sigmaInv = 1/sigma;
	  s.scale(sigmaInv);
	  s.axpy(r, sigmaInv);
	  delta = s.dotRowRow(0, p, 0);
	}
      
      // 3 Scale s_k
      double lambdaDiff = lambda-lambdaBar;
      s.axpy(p, lambdaDiff);
      delta += lambdaDiff*normp;
      
      // 4 
      if(delta <= 0.0) // Make Hessian positive definite.
	{
	  double deltaOverNormp2 = delta/normp2;
	  s.axpy(p, (lambda-2.0*deltaOverNormp2));
	  
	  lambdaBar = 2.0*(lambda - deltaOverNormp2);
	  delta = lambda*normp2 - delta; 
	  lambda = lambdaBar;
	}
      
      // 5  Calculate step size.
      mu=p.dotRowRow(0, r, 0);
      alpha = mu/delta;
      
      // 6 Compute the comparision parameter.
      wPlus.deepCopy(w);
      wPlus.axpy(p, alpha);
      setOptParams(wPlus);
      newObj = computeObjectiveVal();
      Delta = 2.0*delta*(oldObj - newObj)/(mu*mu);

      // 7 Check whether a successful error reduction can be made.
      if(Delta >= 0.0)  // update is successful
	{
	  w.deepCopy(wPlus); 	  
	  oldObj = newObj;
	  computeObjectiveGradParams(rp);
	  rp.negate();
	  lambdaBar = 0; 
	  success = true;
	  // 7.a Check for algorithm restart.
	  if(k % nParams == 0) // restart algorithm
	    p.deepCopy(rp); 
	  else
	    {	      
	      double rpnorm2 = rp.norm2Row(0);
	      double rrp = r.dotRowRow(0, rp, 0); 
	      
	      beta = (rpnorm2 - rrp)/mu;
	      p.scale(beta); 
	      p.axpy(rp, 1.0);
	    }
	  
	  r.deepCopy(rp);
	  
	  // 7.b Reduce the scale parameter
	  if(Delta >= 0.75) lambda *= 0.5;
	  if(lambda<1e-15) lambda = 1e-15;
	}
      else // no reduction in error is possible.
	{
	  setOptParams(w);
	  lambdaBar = lambda; 
	  success = false; 
	}
      
      // 8 Increase the scale parameter
      if(Delta < 0.25) lambda *= 4.0;
      
      // 9 Check for convergence       
      if(getVerbosity()>2)
	cout << "Iteration: " << k << " Error: " << oldObj << " Scale: " << lambda << endl;
      if (success && abs(p.max()*alpha) < paramTol && max(abs(newObj-oldObj)) < objectiveTol)
	return;
      

    }
  cout << "Warning: Maximum number of iterations has been exceeded" << endl;
}
