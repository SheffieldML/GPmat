function k = whitehKernDiagCompute(kern, x)

% WHITEHKERNDIAGCOMPUTE Compute diagonal of WHITEH kernel.
%
%	Description:
%
%	K = WHITEHKERNDIAGCOMPUTE(KERN, X) computes the diagonal of the
%	kernel matrix for the whiteh noise kernel given a design matrix of
%	inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	WHITEHKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, WHITEHKERNCOMPUTE


%	Copyright (c) Neil D. Lawrence, 2009 Mauricio A. Alvarez


k = kern.variance./x;

