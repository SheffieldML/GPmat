function K = ggKernDiagCompute(kern, x)

% GGKERNDIAGCOMPUTE Compute diagonal of GG kernel.
%
%	Description:
%
%	K = GGKERNDIAGCOMPUTE(KERN) computes the diagonal of the kernel
%	matrix for the gaussian gaussian kernel given a design matrix of
%	inputs.
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
%	GGKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, GGKERNCOMPUTE


%	Copyright (c) 2008 Mauricio A. Alvarez and Neil D. Lawrence


%	With modifications by Mauricio A. Alvarez 2009


K = kern.sigma2Latent*kern.sensitivity^2*ones(size(x,1),1);

