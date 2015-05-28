function k = linard2KernDiagCompute(kern, x)

% LINARD2KERNDIAGCOMPUTE Compute diagonal of LINARD2 kernel.
%
%	Description:
%
%	K = LINARD2KERNDIAGCOMPUTE(KERN, X) computes the diagonal of the
%	kernel matrix for the automatic relevance determination linear
%	kernel given a design matrix of inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%	
%
%	See also
%	LINARD2KERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, LINARD2KERNCOMPUTE


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence
%	Copyright (c) 2009 Michalis K. Titsias



scales = sparse(diag(sqrt(kern.inputScales)));
x = x*scales;

k = sum(x.*x, 2);
