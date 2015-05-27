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
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009.
  
% KERN

k = kern.sigma2Noise*ones(size(x,1),1);

