function g = ggKernDiagGradient(kern, x, covDiag)

% GGKERNDIAGGRADIENT Compute gradient of the diagonal of GG kernel.
% FORMAT
% DESC computes the gradient of the diagonal of the kernel matrix for the 
% gaussian gaussian kernel given a design matrix of inputs.
% RETURN g  : a vector containing the gradients
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% ARG covDiag : partial derivative wrt diagonal of the covariance matrix
%	
% SEEALSO : ggKernParamInit, kernDiagCompute, ggKernDiagCompute
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN


grad_sigma2Latent = (kern.sensitivity^2)*sum(covDiag);
grad_sensitivity = (2*kern.sigma2Latent*kern.sensitivity)*sum(covDiag);


g = [zeros(size(kern.precisionU(:)')) zeros(size(kern.precisionG(:)')) ...
    grad_sigma2Latent grad_sensitivity];

