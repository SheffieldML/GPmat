function [k, n2] = rbfard2KernCompute(kern, x, x2)

% RBFARD2KERNCOMPUTE Compute the RBFARD kernel given the parameters and X.
%
%	Description:
%
%	K = RBFARD2KERNCOMPUTE(KERN, X, X2) computes the kernel parameters
%	for the automatic relevance determination radial basis function
%	kernel given inputs associated with rows and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the input matrix associated with the columns of the kernel.
%
%	K = RBFARD2KERNCOMPUTE(KERN, X) computes the kernel matrix for the
%	automatic relevance determination radial basis function kernel given
%	a design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%	
%
%	See also
%	RBFARD2KERNPARAMINIT, KERNCOMPUTE, KERNCREATE, RBFARD2KERNDIAGCOMPUTE


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence
%	Copyright (c) 2009 Michalis K. Titsias



scales = sparse(diag(sqrt(kern.inputScales)));
x = x*scales;
    
if nargin < 3
  n2 = dist2(x, x);
  k = kern.variance*exp(-n2*0.5);
else
  x2 = x2*scales;
  n2 = dist2(x, x2);
  k = kern.variance*exp(-n2*0.5);
end
