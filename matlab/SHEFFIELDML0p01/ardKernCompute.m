function [k, rbfPart, linearPart, n2] = ardKernCompute(kern, x, x2)

% ARDKERNCOMPUTE Compute the ARD kernel given the parameters and X.
%
%	Description:
%
%	K = ARDKERNCOMPUTE(KERN, X, X2) computes the kernel parameters for
%	the pre-built RBF and linear ARD kernel given inputs associated with
%	rows and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the input matrix associated with the columns of the kernel.
%
%	K = ARDKERNCOMPUTE(KERN, X) computes the kernel matrix for the
%	pre-built RBF and linear ARD kernel given a design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	ARDKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, ARDKERNDIAGCOMPUTE


%	Copyright (c) 2004 Neil D. Lawrence



scales = diag(sqrt(kern.inputScales));
x = x*scales;
    
if nargin < 3
  n2 = dist2(x, x);
  wi2 = (.5 .* kern.inverseWidth);
  rbfPart = kern.rbfVariance*exp(-n2*wi2);
  linearPart = x*x'*kern.linearVariance;
  k = rbfPart + kern.whiteVariance*eye(size(x, 1)) + linearPart;
else
  x2 = x2*scales;
  n2 = dist2(x, x2);
  wi2 = (.5 .* kern.inverseWidth);
  rbfPart = kern.rbfVariance*exp(-n2*wi2);
  linearPart = x*x2'*kern.linearVariance;
  k = rbfPart + linearPart;
end
k = k + kern.biasVariance;
