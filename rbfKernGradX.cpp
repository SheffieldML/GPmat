#include "mex.h"
#include <cmath>

void rbfKernGradX(double gX[], 
		  double X[], 
		  double X2[], 
		  int Xrows,
		  int Xcols,
		  int X2rows,
		  double inverseWidth,
		  double variance)

{
  double wi2 = .5*inverseWidth;
  double pf = variance*inverseWidth;
  int ind = 0;
  for(int i = 0; i < X2rows; i++)
    {
      for(int k = 0; k < Xrows; k++)
	{
	  // compute distance between two vectors
	  double n2 = 0;	
	  for(int j = 0; j < Xcols; j++)
	    {	      
	      ind = i+j*X2rows+k*X2rows*Xcols;
	      gX[ind] = X2[i+j*X2rows]-X[k+j*Xrows];
	      n2 += gX[ind]*gX[ind];
	    }
	  // compute gradient of kernel.
	  for(j = 0; j < Xcols; j++)
	    {
	      ind = i+j*X2rows+k*X2rows*Xcols;
	      gX[ind] = pf
		*gX[ind]
		*exp(-n2*wi2);
	    
	    }
	}
    }
}

void mexFunction(
                 int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

  // check proper input and output
  if(nrhs!=3)
    mexErrMsgTxt("Three inputs required.");
  if(nlhs > 1)
    mexErrMsgTxt("Too many output arguments.");
  
  // The kernel structure.
  if(mxGetClassID(prhs[0]) != mxSTRUCT_CLASS)
    mexErrMsgTxt("Error kern should be STRUCT");  
  double* inverseWidth = mxGetPr(mxGetField(prhs[0], 0, "inverseWidth"));
  double* variance = mxGetPr(mxGetField(prhs[0], 0, "variance"));

  // First matrix input.
  if(mxGetClassID(prhs[1]) != mxDOUBLE_CLASS)
    mexErrMsgTxt("Error X should be a numeric array");  
  double* X = mxGetPr(prhs[1]);
  int Xrows = mxGetM(prhs[1]);
  int Xcols = mxGetN(prhs[1]);

  // Second matrix input
  if(mxGetClassID(prhs[2]) != mxDOUBLE_CLASS)
    mexErrMsgTxt("Error X2 should be a numeric array");  
  double* X2 = mxGetPr(prhs[2]);
  int X2rows = mxGetM(prhs[2]);
  int X2cols = mxGetN(prhs[2]);
  if (X2cols != Xcols)
    mexErrMsgTxt("Error number of columns in X and X2 should match");

  // create output.
  int dims[3];
  dims[0] = X2rows;
  dims[1] = Xcols;
  dims[2] = Xrows;
  plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  double* gX = mxGetPr(plhs[0]);

  rbfKernGradX(gX, X, X2, Xrows, Xcols, X2rows, inverseWidth[0], variance[0]);
}
