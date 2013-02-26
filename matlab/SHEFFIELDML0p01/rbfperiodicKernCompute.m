function [k, sk, n2] = rbfperiodicKernCompute(kern, x, x2)

% RBFPERIODICKERNCOMPUTE Compute the RBFPERIODIC kernel given the parameters and X.
%
%	Description:
%
%	K = RBFPERIODICKERNCOMPUTE(KERN, X, X2) computes the kernel
%	parameters for the RBF derived periodic kernel given inputs
%	associated with rows and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the input matrix associated with the columns of the kernel.
%
%	K = RBFPERIODICKERNCOMPUTE(KERN, X) computes the kernel matrix for
%	the RBF derived periodic kernel given a design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	RBFPERIODICKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, RBFPERIODICKERNDIAGCOMPUTE


%	Copyright (c) 2007, 2009 Neil D. Lawrence


factor = 2*pi/kern.period;
if nargin < 3
  n2 = sin(0.5*factor*(repmat(x, 1, size(x, 1)) - repmat(x', size(x, 1), 1)));
  n2 = n2.*n2;
  wi2 = (2 .* kern.inverseWidth);
  sk = exp(-n2*wi2);
else
  n2 = sin(0.5*factor*(repmat(x, 1, size(x2, 1)) - repmat(x2', size(x, 1), 1)));  
  n2 = n2.*n2;
  wi2 = (2 .* kern.inverseWidth);
  sk = exp(-n2*wi2);
end
k = kern.variance*sk;
  