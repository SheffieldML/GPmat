function k = biasKernDiagCompute(kern, x)

% BIASKERNDIAGCOMPUTE Compute diagonal of BIAS kernel.
%
%	Description:
%
%	K = BIASKERNDIAGCOMPUTE(KERN, X) computes the diagonal of the kernel
%	matrix for the bias kernel given a design matrix of inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	BIASKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, BIASKERNCOMPUTE


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



k = repmat(kern.variance, size(x, 1), 1);