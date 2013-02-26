function K = dexpKernDiagCompute(kern, x)

% DEXPKERNDIAGCOMPUTE Compute diagonal of the double exponential kernel.
%
%	Description:
%	
%
%	K = DEXPKERNDIAGCOMPUTE(KERN, X) computes the diagonal of the kernel
%	matrix for the double exponential kernel given a column vector of
%	inputs.
%	 Returns:
%	  K - a vector of the same size as x containing the diagonal of the
%	   kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the kernel matrix is
%	   computed.
%	  X - input data in the form of a design matrix.
%	
%
%	See also
%	DEXPKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, DEXPKERNCOMPUTE


%	Copyright (c) 2009 David Luengo



K = 0.5 * kern.variance * kern.decay * ones(size(x));
