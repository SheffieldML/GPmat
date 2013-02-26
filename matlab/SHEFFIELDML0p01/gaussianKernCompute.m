function [K, Kbase, n2]  = gaussianKernCompute(kern, x, x2)

% GAUSSIANKERNCOMPUTE Compute the Gaussian kernel given the parameters and X.
%
%	Description:
%
%	K = GAUSSIANKERNCOMPUTE(KERN, X, X2) computes the kernel parameters
%	for the Gaussian kernel given inputs associated with rows and
%	columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the input matrix associated with the columns of the kernel.
%
%	K = GAUSSIANKERNCOMPUTE(KERN, X) computes the kernel matrix for the
%	Gaussian kernel given a design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%	
%
%	See also
%	GAUSSIANKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, GAUSSIANKERNDIAGCOMPUTE


%	Copyright (c) 2008 Mauricio Alvarez and Neil D. Lawrence


%	With modifications by Mauricio Alvarez 2009


if kern.isArd
    sqrtP = sqrt(kern.precisionU);
    sqrtPx = x*sparseDiag(sqrtP);
    if nargin < 3
        n2 = dist2(sqrtPx, sqrtPx);        
    else
        sqrtPx2 = x2*sparseDiag(sqrtP);
        n2 = dist2(sqrtPx, sqrtPx2);        
    end
    Kbase = exp(-0.5*n2);    
else
    if nargin < 3
        n2 = dist2(x, x);        
    else        
        n2 = dist2(x, x2);        
    end
    Kbase = exp(-0.5*kern.precisionU*n2);    
end
K = kern.sigma2Latent*Kbase;    








