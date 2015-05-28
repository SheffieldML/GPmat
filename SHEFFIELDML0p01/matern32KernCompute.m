function [K, sk, sqrtn2wi2] = matern32KernCompute(kern, x, x2)

% MATERN32KERNCOMPUTE Compute the MATERN32 kernel given the parameters and X.
%
%	Description:
%
%	K = MATERN32KERNCOMPUTE(KERN, X, X2) computes the kernel parameters
%	for the matern kernel with nu=3/2 kernel given inputs associated
%	with rows and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the input matrix associated with the columns of the kernel.
%
%	K = MATERN32KERNCOMPUTE(KERN, X) computes the kernel matrix for the
%	matern kernel with nu=3/2 kernel given a design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	MATERN32KERNPARAMINIT, KERNCOMPUTE, KERNCREATE, MATERN32KERNDIAGCOMPUTE


%	Copyright (c) 2006, 2009 Neil D. Lawrence


if nargin < 3
  n2 = dist2(x, x);
else
  n2 = dist2(x, x2);
end
wi2 = (3/(kern.lengthScale*kern.lengthScale));
sqrtn2wi2 = sqrt(n2*wi2);
sk = (1+sqrtn2wi2).*exp(-sqrtn2wi2);
K = kern.variance*sk;
