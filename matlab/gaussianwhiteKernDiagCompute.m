function k = gaussianwhiteKernDiagCompute(kern, x)

% GAUSSIANWHITEKERNDIAGCOMPUTE Compute diagonal of gaussian white kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the gaussian white kernel
% given a design matrix of inputs.
% RETURN K : a vector containing the diagonal of the kernel matrix computed
% at the given points.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG X : input data matrix in the form of a design matrix.
%	
% SEEALSO : gaussianwhiteKernParamInit, kernDiagCompute, kernCreate,
% gaussianwhiteKernCompute
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008
  
% KERN
  
Pinv = 2./kern.precisionT;
factor = kern.sigma2Noise/((2*pi)^(kern.inputDimension/2)*sqrt(prod(Pinv))); 

k = factor*ones(size(x,1),1);