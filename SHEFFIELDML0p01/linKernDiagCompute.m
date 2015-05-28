function k = linKernDiagCompute(kern, x)

% LINKERNDIAGCOMPUTE Compute diagonal of LIN kernel.
%
%	Description:
%
%	K = LINKERNDIAGCOMPUTE(KERN, X) computes the diagonal of the kernel
%	matrix for the linear kernel given a design matrix of inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	LINKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, LINKERNCOMPUTE


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



k =  sum(x.*x, 2)*kern.variance;
