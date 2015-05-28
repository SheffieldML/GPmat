function [k, sk] = biasKernCompute(kern, x, x2)

% BIASKERNCOMPUTE Compute the BIAS kernel given the parameters and X.
%
%	Description:
%
%	K = BIASKERNCOMPUTE(KERN, X, X2) computes the kernel parameters for
%	the bias kernel given inputs associated with rows and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the inpute matrix associated with the columns of the kernel.
%
%	K = BIASKERNCOMPUTE(KERN, X) computes the kernel matrix for the bias
%	kernel given a design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	BIASKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, BIASKERNDIAGCOMPUTE


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



if nargin< 3
  dims = [size(x, 1), size(x, 1)];
else
  dims = [size(x, 1), size(x2, 1)];
end
k = repmat(kern.variance, dims);
if nargout > 1
  sk = ones(dims);
end