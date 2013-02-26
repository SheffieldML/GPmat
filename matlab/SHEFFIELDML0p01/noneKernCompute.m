function k = noneKernCompute(kern, x, x2)

% NONEKERNCOMPUTE Compute the NONE kernel given the parameters and X.
%
%	Description:
%
%	K = NONEKERNCOMPUTE(KERN, X, X2) computes the kernel parameters for
%	the dummy kernel function kernel given inputs associated with rows
%	and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the input matrix associated with the columns of the kernel.
%
%	K = NONEKERNCOMPUTE(KERN, X) computes the kernel matrix for the
%	dummy kernel function kernel given a design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	NONEKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, NONEKERNDIAGCOMPUTE


%	Copyright (c) 2008 Neil D. Lawrence




if nargin < 3
  k = spalloc(size(x, 1), size(x, 1), 0);
else
  k = spalloc(size(x, 1), size(x2, 1), 0);
end
