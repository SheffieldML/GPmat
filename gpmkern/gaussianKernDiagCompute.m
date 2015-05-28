function k = gaussianKernDiagCompute(kern, x)

% GAUSSIANKERNDIAGCOMPUTE Compute diagonal of gaussian kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the gaussian kernel
% given a design matrix of inputs.
% RETURN K : a vector containing the diagonal of the kernel matrix computed
% at the given points.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG X : input data matrix in the form of a design matrix.
%	
% SEEALSO : gaussianKernParamInit, kernDiagCompute, kernCreate,
% gaussianKernCompute
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009
  
% KERN
  
k = kern.sigma2Latent*ones(size(x,1),1);
