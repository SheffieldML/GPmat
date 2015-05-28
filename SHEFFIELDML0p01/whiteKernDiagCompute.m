function k = whiteKernDiagCompute(kern, x)

% WHITEKERNDIAGCOMPUTE Compute diagonal of WHITE kernel.
%
%	Description:
%
%	K = WHITEKERNDIAGCOMPUTE(KERN, X) computes the diagonal of the
%	kernel matrix for the white noise kernel given a design matrix of
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
%	WHITEKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, WHITEKERNCOMPUTE


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



k = repmat(kern.variance, size(x, 1), 1);
