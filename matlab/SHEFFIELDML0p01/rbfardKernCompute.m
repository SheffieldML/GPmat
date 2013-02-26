function [k, n2] = rbfardKernCompute(kern, x, x2)

% RBFARDKERNCOMPUTE Compute the RBFARD kernel given the parameters and X.
%
%	Description:
%
%	K = RBFARDKERNCOMPUTE(KERN, X, X2) computes the kernel parameters
%	for the automatic relevance determination radial basis function
%	kernel given inputs associated with rows and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the input matrix associated with the columns of the kernel.
%
%	K = RBFARDKERNCOMPUTE(KERN, X) computes the kernel matrix for the
%	automatic relevance determination radial basis function kernel given
%	a design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	RBFARDKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, RBFARDKERNDIAGCOMPUTE


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



scales = sparse(diag(sqrt(kern.inputScales)));
x = x*scales;
    
if nargin < 3
  n2 = dist2(x, x);
  wi2 = (.5 .* kern.inverseWidth);
  k = kern.variance*exp(-n2*wi2);
else
  x2 = x2*scales;
  n2 = dist2(x, x2);
  wi2 = (.5 .* kern.inverseWidth);
  k = kern.variance*exp(-n2*wi2);
end
