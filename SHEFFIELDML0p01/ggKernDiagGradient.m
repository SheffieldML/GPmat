function g = ggKernDiagGradient(kern, x, covDiag)

% GGKERNDIAGGRADIENT Compute gradient of the diagonal of GG kernel.
%
%	Description:
%
%	G = GGKERNDIAGGRADIENT(KERN, X, COVDIAG) computes the gradient of
%	the diagonal of the kernel matrix for the gaussian gaussian kernel
%	given a design matrix of inputs.
%	 Returns:
%	  G - a vector containing the gradients
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	  COVDIAG - partial derivative wrt diagonal of the covariance matrix
%	
%
%	See also
%	GGKERNPARAMINIT, KERNDIAGCOMPUTE, GGKERNDIAGCOMPUTE


%	Copyright (c) 2010 Mauricio A. Alvarez



grad_sigma2Latent = (kern.sensitivity^2)*sum(covDiag);
grad_sensitivity = (2*kern.sigma2Latent*kern.sensitivity)*sum(covDiag);


g = [zeros(size(kern.precisionU(:)')) zeros(size(kern.precisionG(:)')) ...
    grad_sigma2Latent grad_sensitivity];

