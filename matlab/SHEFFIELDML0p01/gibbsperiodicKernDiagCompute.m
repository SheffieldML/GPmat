function k = gibbsperiodicKernDiagCompute(kern, x)

% GIBBSPERIODICKERNDIAGCOMPUTE Compute diagonal of GIBBSPERIODIC kernel.
%
%	Description:
%
%	K = GIBBSPERIODICKERNDIAGCOMPUTE(KERN, X) computes the diagonal of
%	the kernel matrix for the Gibbs-kernel derived periodic kernel given
%	a design matrix of inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	GIBBSPERIODICKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, GIBBSPERIODICKERNCOMPUTE


%	Copyright (c) 2007 Neil D. Lawrence


k = repmat(kern.variance, size(x, 1), 1);

