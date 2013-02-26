function k = rbfperiodic2KernDiagCompute(kern, x)

% RBFPERIODIC2KERNDIAGCOMPUTE Compute diagonal of RBFPERIODIC2 kernel.
%
%	Description:
%
%	K = RBFPERIODIC2KERNDIAGCOMPUTE(KERN, X) computes the diagonal of
%	the kernel matrix for the RBF periodic covariance with variying
%	period kernel given a design matrix of inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%	
%	
%
%	See also
%	RBFPERIODIC2KERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, RBFPERIODIC2KERNCOMPUTE


%	Copyright (c) 2007, 2009 Neil D. Lawrence


%	With modifications by Andreas C. Damianou 2011


%	With modifications by Michalis K. Titsias 2011




k = repmat(kern.variance, size(x, 1), 1);
