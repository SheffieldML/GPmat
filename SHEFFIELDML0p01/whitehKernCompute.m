function [k, sk] = whitehKernCompute(kern, x, x2)

% WHITEHKERNCOMPUTE Compute the WHITEH kernel given the parameters and X.
%
%	Description:
%
%	K = WHITEHKERNCOMPUTE(KERN, X, X2) computes the kernel parameters
%	for the whiteh noise kernel given inputs associated with rows and
%	columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the inpute matrix associated with the columns of the kernel.
%
%	K = WHITEHKERNCOMPUTE(KERN, X) computes the kernel matrix for the
%	whiteh noise kernel given a design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	WHITEHKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, WHITEHKERNDIAGCOMPUTE


%	Copyright (c) Neil D. Lawrence, 2009 Mauricio A. Alvarez


if nargin < 3
  sk = sparseDiag(1./x(:,end));  
  k = kern.variance*sk;
else
  k = spalloc(size(x, 1), size(x2, 1), 0);
end
