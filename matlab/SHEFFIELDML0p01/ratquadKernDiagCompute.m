function k = ratquadKernDiagCompute(kern, x)

% RATQUADKERNDIAGCOMPUTE Compute diagonal of RATQUAD kernel.
%
%	Description:
%
%	K = RATQUADKERNDIAGCOMPUTE(KERN, X) computes the diagonal of the
%	kernel matrix for the rational quadratic kernel given a design
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
%	RATQUADKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, RATQUADKERNCOMPUTE


%	Copyright (c) 2006 Neil D. Lawrence


k = repmat(kern.variance, size(x, 1), 1);
