#ifndef COPTIMISABLE_H
#define COPTIMISABLE_H
#include "CMatrix.h"

// abstract class for making a class optimisable.
class COptimisable {

 public:
  virtual int getOptNumParams() const=0;
  virtual void getOptParams(CMatrix& param) const=0;
  virtual void setOptParams(const CMatrix& param)=0;
  virtual void computeObjectiveGradParams(CMatrix& g)=0;
  virtual double computeObjectiveVal()=0;
  void checkGradients()
    {
      double changeFactor = 1e-6;
      double origParam = 0.0;
      double change = 0.0;
      double objectivePlus = 0.0;
      double objectiveMinus = 0.0;

      CMatrix analyticGrad(1, getOptNumParams());
      CMatrix numericalDiff(1, getOptNumParams());
      CMatrix diffNumericalAnalytic(1, getOptNumParams());
      CMatrix params(1, getOptNumParams());
      CMatrix origParams(1, getOptNumParams());

      getOptParams(params);
      origParams.deepCopy(params);
      for(int j=0; j<getOptNumParams(); j++)
	{
	  origParam = origParams.getVals(j);
	  change = changeFactor*origParam;
	  params.setVals(origParam + change, j);
	  setOptParams(params);
	  objectivePlus = computeObjectiveVal();
	  params.setVals(origParam - change, j);
	  setOptParams(params);
	  objectiveMinus = computeObjectiveVal();
	  numericalDiff.setVals(0.5*(objectivePlus - objectiveMinus)/change, j);
	  params.setVals(origParams.getVals(j), j);
	}
  
      cout << "Numerical differences:" << endl << numericalDiff << endl;
      computeObjectiveGradParams(analyticGrad);
      cout << "Analytic gradients:" << endl << analyticGrad << endl;
      diffNumericalAnalytic.deepCopy(analyticGrad);
      diffNumericalAnalytic-=numericalDiff;
      cout << "Maximum param difference: " << diffNumericalAnalytic.max() << endl;
    }
  void gdOptimise(double learnRate, double momentum, int display, int maxIters)
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
	  if(display)
	    {
	      cout << "Iteration: " << iter << ", objective function: " << objectiveVal << endl;
	    }
	  diffParam = params.maxAbsDiff(oldParams);
	  if(diffObjective<objectiveTol && diffParam<parameterTol)
	    {
	      cout << "Converged .." << endl;
	      break;
	    }

	} 
      cout << "Parameters: " << endl;
      cout << params << endl;
/*       if(iter>=maxIters) */
/* 	cout << "Maximum iterations exceed in gdOptimise(), objective change, " << diffObjective << ", parameter change, " << diffParam << endl; */
   }

  
 private:
  
  double objectiveTol;
  double parameterTol;
  //int maxIters;
  //double learnRate;
  //double momentum;
  bool evalFunc;
  //int display;


}; 

#endif
