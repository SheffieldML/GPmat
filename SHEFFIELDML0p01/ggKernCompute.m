function [K, Kbase, Prinv, Pqrinv, P, n2] = ggKernCompute(kern, x, x2)

% GGKERNCOMPUTE Compute the GG kernel given the parameters and X.
%
%	Description:
%
%	K = GGKERNCOMPUTE(KERN, X, KERNEL.) computes the kernel parameters
%	for the gaussian gaussian kernel given inputs associated with rows
%	and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  KERNEL. - % ARG x2 the input matrix associated with the columns of
%	   the kernel.
%
%	K = GGKERNCOMPUTE(KERN, X) computes the kernel matrix for the
%	gaussian kernel given a design matrix of inputs. Returns:
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%	
%
%	See also
%	GGKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, GGKERNDIAGCOMPUTE


%	Copyright (c) 2008 Mauricio A. Alvarez and Neil D. Lawrence


%	With modifications by Mauricio A. Alvarez 2009.



if nargin < 3
  x2 = x;
end

Pr  = kern.precisionU;
Pqr = kern.precisionG;
Prinv = 1./Pr;
Pqrinv = 1./Pqr;
Pinv = Prinv + 2*Pqrinv;
P = 1./Pinv;
if kern.isArd,
    sqrtP = sparseDiag(sqrt(P));
    sqrtPx = x*sqrtP;
    sqrtPx2 = x2*sqrtP;
    n2 = dist2(sqrtPx, sqrtPx2);
    Kbase = exp(-0.5*n2);    
else
    n2 = dist2(x, x2);
    Kbase = exp(-0.5*P*n2);    
end
K = kern.sigma2Latent*kern.sensitivity^2*...
    Kbase;

