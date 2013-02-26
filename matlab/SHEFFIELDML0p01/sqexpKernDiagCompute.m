function k = sqexpKernDiagCompute(kern, x)

% SQEXPKERNDIAGCOMPUTE Compute diagonal of SQEXP kernel.
%
%	Description:
%
%	K = SQEXPKERNDIAGCOMPUTE(KERN, X) computes the diagonal of the
%	kernel matrix for the pre-built compound squared exponential kernel
%	given a design matrix of inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	SQEXPKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, SQEXPKERNCOMPUTE


%	Copyright (c) 2004 Neil D. Lawrence



k = repmat(kern.rbfVariance+kern.whiteVariance, size(x, 1), 1)  + kern.biasVariance;
