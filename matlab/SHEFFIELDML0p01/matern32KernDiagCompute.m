function k = matern32KernDiagCompute(kern, x)

% MATERN32KERNDIAGCOMPUTE Compute diagonal of MATERN32 kernel.
%
%	Description:
%
%	K = MATERN32KERNDIAGCOMPUTE(KERN, X) computes the diagonal of the
%	kernel matrix for the matern kernel with nu=3/2 kernel given a
%	design matrix of inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	MATERN32KERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, MATERN32KERNCOMPUTE


%	Copyright (c) 2006 Neil D. Lawrence


k = repmat(kern.variance, size(x, 1), 1);
