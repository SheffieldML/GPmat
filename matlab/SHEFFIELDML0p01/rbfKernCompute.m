function [k, sk, n2] = rbfKernCompute(kern, x, x2)

% RBFKERNCOMPUTE Compute the RBF kernel given the parameters and X.
%
%	Description:
%
%	[K, SK] = RBFKERNCOMPUTE(KERN, X, X2) computes the kernel parameters
%	for the radial basis function kernel given inputs associated with
%	rows and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	  SK - unscaled kernel matrix (i.e. only the exponential part).
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the input matrix associated with the columns of the kernel.
%
%	[K, SK] = RBFKERNCOMPUTE(KERN, X) computes the kernel matrix for the
%	radial basis function kernel given a design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	  SK - unscaled kernel matrix (i.e. only the exponential part).
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%	
%
%	See also
%	RBFKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, RBFKERNDIAGCOMPUTE


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence


%	With modifications by Mauricio Alvarez 2009, David Luengo, 2009


if nargin < 3
  n2 = dist2(x, x);
else
  n2 = dist2(x, x2);
end

wi2 = (.5 .* kern.inverseWidth);
sk = exp(-n2*wi2);
k = kern.variance*sk;
if isfield(kern, 'isNormalised') && (kern.isNormalised == true)
    k = k * sqrt(kern.inverseWidth/(2*pi));
end
