function k = ggwhiteKernDiagCompute(kern, x)

% GGWHITEKERNDIAGCOMPUTE Compute diagonal of GG WHITE kernel.
%
%	Description:
%
%	K = GGWHITEKERNDIAGCOMPUTE(KERN) computes the diagonal of the kernel
%	matrix for the gaussian gaussian white kernel given a design matrix
%	of inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed. X :
%	   input data matrix in the form of a design matrix.
%	
%	
%
%	See also
%	GGWHITEKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, GGWHITEKERNCOMPUTE


%	Copyright (c) 2008 Mauricio A. Alvarez and Neil D. Lawrence


%	With modifications by Mauricio A. Alvarez 2009.


k = kern.sigma2Noise*kern.variance^2*ones(size(x,1),1);