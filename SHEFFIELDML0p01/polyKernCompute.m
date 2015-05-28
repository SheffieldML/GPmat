function k = polyKernCompute(kern, x, x2)

% POLYKERNCOMPUTE Compute the POLY kernel given the parameters and X.
%
%	Description:
%
%	K = POLYKERNCOMPUTE(KERN, X, X2) computes the kernel parameters for
%	the polynomial kernel given inputs associated with rows and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the input matrix associated with the columns of the kernel.
%
%	K = POLYKERNCOMPUTE(KERN, X) computes the kernel matrix for the
%	polynomial kernel given a design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	POLYKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, POLYKERNDIAGCOMPUTE


%	Copyright (c) 2005, 2006 Neil D. Lawrence



if nargin < 3
  k = kern.variance*(x*x'*kern.weightVariance+kern.biasVariance).^kern.degree;
else
  k = kern.variance*(kern.weightVariance*x*x2'+kern.biasVariance).^kern.degree;
end
if issparse(x)
  k = full(k);
end