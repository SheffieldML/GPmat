function [k, sk] = whitefixedKernCompute(kern, x, x2)

% WHITEFIXEDKERNCOMPUTE Compute the WHITEFIXED kernel given the parameters and X.
%
%	Description:
%
%	K = WHITEFIXEDKERNCOMPUTE(KERN, X, X2) computes the kernel
%	parameters for the fixed parameter white noise kernel given inputs
%	associated with rows and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the inpute matrix associated with the columns of the kernel.
%
%	K = WHITEFIXEDKERNCOMPUTE(KERN, X) computes the kernel matrix for
%	the fixed parameter white noise kernel given a design matrix of
%	inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%	
%
%	See also
%	WHITEFIXEDKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, WHITEFIXEDKERNDIAGCOMPUTE


%	Copyright (c) 2006 Nathaniel J. King


%	With modifications by Neil D. Lawrence 2009



if nargin < 3
  k = whiteKernCompute(kern, x);
else
  k = whiteKernCompute(kern, x, x2);
end
if nargout > 1 
  sk = k;
end