function [k, sk] = linKernCompute(kern, x, x2)

% LINKERNCOMPUTE Compute the LIN kernel given the parameters and X.
%
%	Description:
%
%	K = LINKERNCOMPUTE(KERN, X, X2) computes the kernel parameters for
%	the linear kernel given inputs associated with rows and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the input matrix associated with the columns of the kernel.
%
%	K = LINKERNCOMPUTE(KERN, X) computes the kernel matrix for the
%	linear kernel given a design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	LINKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, LINKERNDIAGCOMPUTE


%	Copyright (c) 2004, 2005, 2006, 2009 Neil D. Lawrence



if nargin < 3
  sk = x*x'; 
else
  sk = x*x2';
end
k = sk*kern.variance;
if issparse(x)
  k = full(k);
end
