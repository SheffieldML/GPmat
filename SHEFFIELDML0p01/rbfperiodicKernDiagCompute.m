function k = rbfperiodicKernDiagCompute(kern, x)

% RBFPERIODICKERNDIAGCOMPUTE Compute diagonal of RBFPERIODIC kernel.
%
%	Description:
%
%	K = RBFPERIODICKERNDIAGCOMPUTE(KERN, X) computes the diagonal of the
%	kernel matrix for the RBF derived periodic kernel given a design
%	matrix of inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	RBFPERIODICKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, RBFPERIODICKERNCOMPUTE


%	Copyright (c) 2007 Neil D. Lawrence



k = repmat(kern.variance, size(x, 1), 1);
