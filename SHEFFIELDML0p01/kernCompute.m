function k = kernCompute(kern, x, x2)

% KERNCOMPUTE Compute the kernel given the parameters and X.
%
%	Description:
%
%	K = KERNCOMPUTE(KERN, X) computes a kernel matrix for the given
%	kernel type given an input data matrix.
%	 Returns:
%	  K - computed elements of the kernel structure.
%	 Arguments:
%	  KERN - kernel structure to be computed.
%	  X - input data matrix (rows are data points) to the kernel
%	   computation.
%
%	K = KERNCOMPUTE(KERN, X, X2) computes a kernel matrix for the given
%	kernel type given two input data matrices, one for the rows and one
%	for the columns.
%	 Returns:
%	  K - computed elements of the kernel structure.
%	 Arguments:
%	  KERN - kernel structure to be computed.
%	  X - first input matrix to the kernel computation (forms the rows
%	   of the kernel).
%	  X2 - second input matrix to the kernel computation (forms the
%	   columns of the kernel).
%	
%
%	See also
%	KERNCREATE, KERNDIAGCOMPUTE


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence


fhandle = str2func([kern.type 'KernCompute']);
if nargin < 3
  k = fhandle(kern, x);
else
  k = fhandle(kern, x, x2);
end
