function [k, innerProd, arg, denom, numer, vec] = polyardKernCompute(kern, x, x2)

% POLYARDKERNCOMPUTE Compute the POLYARD kernel given the parameters and X.
%
%	Description:
%
%	K = POLYARDKERNCOMPUTE(KERN, X, X2) computes the kernel parameters
%	for the automatic relevance determination polynomial kernel given
%	inputs associated with rows and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the input matrix associated with the columns of the kernel.
%
%	K = POLYARDKERNCOMPUTE(KERN, X) computes the kernel matrix for the
%	automatic relevance determination polynomial kernel given a design
%	matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	POLYARDKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, POLYARDKERNDIAGCOMPUTE


%	Copyright (c) 2005, 2006 Neil D. Lawrence



scales = sparse(diag(sqrt(kern.inputScales)));
x = x*scales;

if nargin < 3
  innerProd = x*x';
  arg = innerProd*kern.weightVariance + kern.biasVariance;
  k = kern.variance*arg.^kern.degree;
else
  x2 = x2*scales;
  innerProd = x*x2';  
  arg = innerProd*kern.weightVariance + kern.biasVariance;
  k = kern.variance*arg.^kern.degree;
end
