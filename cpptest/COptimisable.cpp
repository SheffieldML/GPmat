#include "COptimisable.h"

void COptimisable::checkGradients()
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
      origParam = origParams.getVal(j);
      change = changeFactor*origParam;
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
  computeObjectiveGradParams(analyticGrad);
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
void COptimisable::bracketMinimum(double& a, double& b, double& c, double& fa, const int maxStep)
{
  // copyright Ian T. Nabney (1996-2001)
  bool terminate=false;
  double fb=oneDObjectiveVal(b);
  double fu;
  double r;
  double q;
  double u;
  double ulimit;
  double fc;
  double denom;
  if(fb>fa)
    {
      do
	{
	  // minimum is between a and b.
	  b = a+(b-a)/phi;
	  c = b;
	  fb=oneDObjectiveVal(b);
	} while(fb>fa);
    }
  else
    {
      c = b+(b-a)*phi;
      fc = oneDObjectiveVal(c);
      while(fb>fc)
	{
	  r = (fb-fc)*(b-a);
	  q = (fb-fa)*(b-c);
	  denom = q-r;
	  if (abs(q-r)<0)
	    denom = smallNum;
	  u = b - ((b-c)*q - (b-a)*r)/(2.0*denom);
	  ulimit = b + maxStep*(c-b);
	  if((b-u)*(u-c)>0.0)
	    {
	      // it is between b and c
	      fu = oneDObjectiveVal(u);
	      if(fu<fc)
		{
		  // minimum between b and c.
		  a = b;
		  b = u;
		  return;
		}
	      else if(fu>fb)
		{
		  // minimum between a and u
		  b = c;
		  c = u;
		  return;
		}		
	      // interpolation with parabola failed, use golden section
	      u = c + (c-b)*phi;
	    }
	  else if((c - u)*(u-ulimit)>0.0)
	    {
	      // it lies between c and ulimit
	      fu = oneDObjectiveVal(u);
	      if(fu<fc)
		{
		  b = c;
		  c = u;
		  u = c + (c-b)*phi;
		}
	      else
		terminate=true;
	    }
	  else if ((u-ulimit)*(ulimit-c)>=0.0)
	    {
	      // limit parabolic u to maximum valie
	      u=ulimit;
	    }
	  else
	    {
	      // reject parabolic u and use golden section step
	      u = c + (c-b)*phi;
	    }
	  if(!terminate)
	    {
	      fu = oneDObjectiveVal(u);
	    }
	  a = b; 
	  b = c;
	  c = u;
	  fa = fb;
	  fb = fc;
	  fc = fu;
	}
    }
  if(a>c)
    {
      double temp = c;
      c = a;
      a = temp;
    } 
}

void COptimisable::lineMinimisation()//double fpt, optionsStuct options)
{
//   // copyright Ian T. Nabney (1996-2001)
//   TOL = sqrt(eps);	// Maximal fractional precision
//   TINY = 1.0e-10;         // Can't use fractional precision when minimum is at 0
 
//  // Bracket the minimum
//   double br_min=0.0;
//   double br_mid=1.0;
//   double br_max=0.0;

//   bracketMinimum(br_min, br_mid, br_max, fpt, const int maxStep)
    
//   // Use Brent's algorithm to find minimum
//   // Initialise the points and function values
//   double w = br_mid;   	// Where second from minimum is
//   double v = br_mid;   	// Previous value of w
//   double x = v;   	// Where current minimum is
//   double e = 0.0; 	// Distance moved on step before last
//   double fx = oneDObjectiveVal(x);
//   double fv = fx; 
//   double fw = fx;
//   double xm = 0.0;
//   double tol1 = 0.0;

//  for(int n=0; n<options.maxIters; n++)
//    {
//      xm = 0.5.*(br_min+br_max);  // Middle of bracket
//      // Make sure that tolerance is big enough
//      tol1 = TOL * (max(abs(x))) + TINY;

//      // check for termination through point position change.
//      if (max(abs(x - xm)) <= options.paramPrecision 
// 	 && br_max-br_min < 4*options.paramPrecision)
//        return;
     
//      // Check if step before last was big enough to try a parabolic step.
//      // Note that this will fail on first iteration, which must be a golden
//      // section step.
//      if (max(abs(e)) > tol1)
//        {
// 	 // Construct a trial parabolic fit through x, v and w
// 	 r = (fx - fv) * (x - w);
// 	 q = (fx - fw) * (x - v);
// 	 p = (x - v)*q - (x - w)*r;
// 	 q = 2.0 * (q - r);
// 	 if (q > 0.0) 
// 	   p = -p;
// 	 q = abs(q);
// 	 // Test if the parabolic fit is OK
// 	 if (abs(p) >= abs(0.5*q*e) || p <= q*(br_min-x) || p >= q*(br_max-x))
// 	   {
// 	     // No it isn't, so take a golden section step
// 	     if (x >= xm)
// 	       e = br_min-x;
// 	     else
// 	       e = br_max-x;
// 	     d = cphi*e;
// 	   }
// 	 else
// 	   {
// 	     // Yes it is, so take the parabolic step
// 	     e = d;
// 	     d = p/q;
// 	     u = x+d;
// 	     if (u-br_min < 2*tol1 || br_max-u < 2*tol1)
// 	       d = sign(xm-x)*tol1;
// 	   }   
//        }
//      else
//        {
// 	 // Step before last not big enough, so take a golden section step
// 	 if (x >= xm)
// 	   e = br_min - x;
// 	 else
// 	   e = br_max - x;
// 	 d = cphi*e;
//        }
//      // Make sure that step is big enough
//      if (abs(d) >= tol1)
//        u = x+d;
//      else
//        u = x + sign(d)*tol1;
//      // Evaluate function at u
//      fu = oneDObjectiveVal(u);
//      // Reorganise bracket
//      if (fu <= fx)
//        {
// 	 if (u >= x)
// 	   br_min = x;
// 	 else
// 	   br_max = x;
// 	 v = w; w = x; x = u;
// 	 fv = fw; fw = fv; fx = fu;
//        }
//      else
//        {
// 	 if (u < x)
// 	   br_min = u;   
// 	 else
// 	   br_max = u;
// 	 if (fu <= fw || w == x)
// 	   {
// 	     v = w; w = u;
// 	     fv = fw; fw = fu;
// 	   }
// 	 else if (fu <= fv || v == x || v == w)
// 	   {
// 	     v = u;
// 	     fv = fu;
// 	   }
//        }
//      if (options.getVerbosity()>2)
//        cout << "Line Minimisation, Iteration: " << n << " Objective function: " << fx << endl;
//    }
//  options(8) = fx;
 
}
void COptimisable::netlabScgOptimise(int maxIters, const double objectiveTol, const double paramTol)
{
  int nParams = getOptNumParams();
  double objectiveVal = 0.0;
  double oldObjective = 0.0;
  double diffObjective = 0.0;
  double diffParam = 0.0;
  
  CMatrix params(1, getOptNumParams());
  CMatrix paramsPlus(1, getOptNumParams());
  CMatrix changeParams(1, getOptNumParams());
  CMatrix paramsNew(1, getOptNumParams());

  CMatrix gradNew(1, getOptNumParams());
  CMatrix gradOld(1, getOptNumParams());
  CMatrix gradChange(1, getOptNumParams());
  CMatrix gradPlus(1, getOptNumParams());
  CMatrix d(1, getOptNumParams());
  changeParams.zeros();
  getOptParams(params);
  if(evalFunc)
    objectiveVal = computeObjectiveVal();

  double sigma0 = 1.0e-4;
  double fold = computeObjectiveVal();
  double fnow = fold;
  computeObjectiveGradParams(gradNew);
  gradOld.deepCopy(gradNew);
  d.deepCopy(gradNew);				// Initial search direction.
  d.negate();
  bool success = true;				// Force calculation of directional derivs.
  int nsuccess = 0;				// nsuccess counts number of successes.
  double beta = 1.0;				// Initial scale parameter.
  double betamin = 1.0e-15; 			// Lower bound on scale.
  double betamax = 1.0e100;			// Upper bound on scale.
  double mu = 0.0;
  double kappa = 0.0;
  double theta = 0.0;
  double sigma = 0.0;
  double delta = 0.0;
  double alpha = 0.0;
  double Delta = 0.0;
  double fnew = 0.0;
  double gamma = 0.0;
  int j = 1;					// j counts number of iterations.
  
  // Main optimization loop.
  while (j <= maxIters)
    {

      // Calculate first and second directional derivatives.
      if (success == 1)
	{
	  for(int i =0; i<d.getRows(); i++)
	    {
	      if(!finite(d.getVal(0, i)))
		{
		  cout << "Warning d is infinite." << endl;
		}    
	    }  
	  mu = d.dotRowRow(0, gradNew, 0);
	  if (mu >= 0)
	    {
	      d.deepCopy(gradNew);
	      d.negate();
	      mu = d.dotRowRow(0, gradNew, 0);
	    }
	  kappa = d.norm2Row(0);
	  if(kappa < ndlutil::EPS)
	    return;
	  sigma = sigma0/sqrt(kappa);
	  paramsPlus.deepCopy(params);
	  paramsPlus.axpy(d, sigma);
	  setOptParams(paramsPlus);
	  computeObjectiveGradParams(gradChange);
	  gradChange.axpy(gradNew, -1.0);
	  theta = gradChange.dotRowRow(0, d, 0)/sigma;
	}
      // Increase effective curvature and evaluate step size alpha.
      delta = theta + beta*kappa;
      if(delta<=0)
	{
	  delta = beta*kappa;
	  beta = beta-theta/kappa;
	}
      alpha=-mu/delta;
      
      // Calculate the comparison ratio.
      paramsNew.deepCopy(params);
      paramsNew.axpy(d, alpha);
      setOptParams(paramsNew);
      fnew = computeObjectiveVal();
      Delta = 2*(fnew - fold)/(alpha*mu);
      if (Delta  >= 0)
	{
	  success = true;
	  nsuccess = nsuccess + 1;
	  params.deepCopy(paramsNew);
	  fnow = fnew;
	}
      else
	{
	  success = false;
	  fnow = fold;
	}
      if(getVerbosity()>2)
	cout << "Iteration: " << j << " Error: " << fnow << " Scale: " << beta << endl;
      
      if (success)
	{
	  // Test for termination
	  if (abs(d.max()*alpha) < paramTol && max(abs(fnew-fold)) < objectiveTol)
	    return;
	
	  else
	    {
	      // Update variables for new position
	      fold = fnew;
	      gradOld.deepCopy(gradNew);
	      setOptParams(params);
	      computeObjectiveGradParams(gradNew);
	      // If the gradient is zero then we are done.
	      if (gradNew.norm2Row(0) == 0)
		return;
	    }
	}

      // Adjust beta according to comparison ratio.
      if (Delta < 0.25)
	beta = min(4.0*beta, betamax);
      if (Delta > 0.75)
	beta = max(0.5*beta, betamin);
      

      // Update search direction using Polak-Ribiere formula, or re-start 
      // in direction of negative gradient after nparams steps.
      if (nsuccess == getOptNumParams())
	{
	  d.deepCopy(gradNew);
	  d.negate();
	  nsuccess = 0;
	}
      else
	if (success)
	  {
	    gradChange.deepCopy(gradOld);
	    gradChange.axpy(gradNew, -1.0);
	    gamma=gradChange.dotRowRow(0, gradNew, 0)/mu;
            if (!finite(gamma))
	      cerr << "gamma is infinity" << endl;
            d.scale(gamma);
	    d.axpy(gradNew, -1.0);
	  }
      
      j++;
    }

  // If we get here, then we haven't terminated in the given number of 
  // iterations.
  cout << "Warning: Maximum number of iterations has been exceeded" << endl;

}
/*	float* w, 
	const ErrFunc* E, 
	unsigned int nits, float errTH)
{
  unsigned int nParams = E->dim();
*/
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
      rr = r.normRow(0); 
      
      if(getVerbosity()>2)
	cout << "Iteration: " << k << " Error: " << newObj << " Scale: " << lambda << endl;
      
      if(rr < paramTol) return;
    }
  cout << "Warning: Maximum number of iterations has been exceeded" << endl;
}
