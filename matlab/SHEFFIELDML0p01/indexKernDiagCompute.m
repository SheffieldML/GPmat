function k = indexKernDiagCompute(kern, x)

% INDEXKERNDIAGCOMPUTE Compute diagonal of INDEX kernel.
%
%	Description:
%
%	K = INDEXKERNDIAGCOMPUTE(KERN, X) computes the diagonal of the
%	kernel matrix for the index based covariance function kernel given a
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
%	INDEXKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, INDEXKERNCOMPUTE


%	Copyright (c) 2011 Neil D. Lawrence

  k = repmat(kern.variance, size(x, 1), 1);
end