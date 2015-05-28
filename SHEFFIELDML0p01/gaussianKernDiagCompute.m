function k = gaussianKernDiagCompute(kern, x)

% GAUSSIANKERNDIAGCOMPUTE Compute diagonal of gaussian kernel.
%
%	Description:
%
%	K = GAUSSIANKERNDIAGCOMPUTE(KERN, X) computes the diagonal of the
%	kernel matrix for the gaussian kernel given a design matrix of
%	inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	gaussianKernCompute
%	
%	
%
%	See also
%	GAUSSIANKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, 


%	Copyright (c) 2008 Mauricio A. Alvarez and Neil D. Lawrence


%	With modifications by Mauricio A. Alvarez 2009

  
k = kern.sigma2Latent*ones(size(x,1),1);
