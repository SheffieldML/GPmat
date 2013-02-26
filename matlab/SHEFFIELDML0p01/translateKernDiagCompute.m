function k = translateKernDiagCompute(kern, x)

% TRANSLATEKERNDIAGCOMPUTE Compute diagonal of TRANSLATE kernel.
%
%	Description:
%
%	K = TRANSLATEKERNDIAGCOMPUTE(KERN, X) computes the diagonal of the
%	kernel matrix for the input space translation kernel given a design
%	matrix of inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	cmpndKernDiagCompute, translateKernCompute
%	
%
%	See also
%	TRANSLATEKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, 


%	Copyright (c) 2007 Neil D. Lawrence


x = x - repmat(kern.centre, size(x, 1), 1);
k = cmpndKernDiagCompute(kern, x);