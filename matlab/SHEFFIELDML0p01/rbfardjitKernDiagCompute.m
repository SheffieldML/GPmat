function k = rbfardjitKernDiagCompute(kern, x)

% RBFARDJITKERNDIAGCOMPUTE Compute diagonal of RBFARDJIT kernel.
%
%	Description:
%
%	K = RBFARDJITKERNDIAGCOMPUTE(KERN, X) computes the diagonal of the
%	kernel matrix for the automatic relevance determination radial basis
%	function kernel given a design matrix of inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%	
%
%	See also
%	RBFARD2KERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, RBFARD2KERNCOMPUTE


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence
%	Copyright (c) 2009 Michalis K. Titsias



k = repmat(kern.variance*(1 + kern.jitter), size(x, 1), 1);
