function k = polyKernDiagCompute(kern, x)

% POLYKERNDIAGCOMPUTE Compute diagonal of POLY kernel.
%
%	Description:
%
%	K = POLYKERNDIAGCOMPUTE(KERN, X) computes the diagonal of the kernel
%	matrix for the polynomial kernel given a design matrix of inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	POLYKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, POLYKERNCOMPUTE


%	Copyright (c) 2005, 2006 Neil D. Lawrence



k =  kern.variance*(sum(x.*x, 2)*kern.weightVariance + kern.biasVariance).^kern.degree;
