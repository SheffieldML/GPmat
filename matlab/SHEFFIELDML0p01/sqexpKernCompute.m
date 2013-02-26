function [k, rbfPart, n2] = sqexpKernCompute(kern, x, x2)

% SQEXPKERNCOMPUTE Compute the SQEXP kernel given the parameters and X.
%
%	Description:
%
%	K = SQEXPKERNCOMPUTE(KERN, X, X2) computes the kernel parameters for
%	the pre-built compound squared exponential kernel given inputs
%	associated with rows and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the input matrix associated with the columns of the kernel.
%
%	K = SQEXPKERNCOMPUTE(KERN, X) computes the kernel matrix for the
%	pre-built compound squared exponential kernel given a design matrix
%	of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	SQEXPKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, SQEXPKERNDIAGCOMPUTE


%	Copyright (c) 2004 Neil D. Lawrence



if nargin < 3
  n2 = dist2(x, x);
  wi2 = (.5 .* kern.inverseWidth);
  rbfPart = kern.rbfVariance*exp(-n2*wi2);
  k = rbfPart + kern.whiteVariance*eye(size(x, 1));
else
  n2 = dist2(x, x2);
  wi2 = (.5 .* kern.inverseWidth);
  rbfPart = kern.rbfVariance*exp(-n2*wi2);
  k = rbfPart;
end
k = k + kern.biasVariance;