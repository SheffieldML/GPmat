function k = ggwhiteKernDiagCompute(kern, x)

% GGWHITEKERNDIAGCOMPUTE Compute diagonal of GG WHITE kernel.
% FORMAT
% DESC computes the diagonal of the kernel
%	matrix for the gaussian gaussian white kernel given a design matrix of
%	inputs.
% RETURN  K : a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
% ARG kern : the kernel structure for which the matrix is computed.
% X : input data matrix in the form of a design matrix.
%	
% SEEALSO : ggwhiteKernParamInit, kernDiagCompute, kernCreate, ggwhiteKernCompute
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009.

% KERN

k = kern.sigma2Noise*kern.variance^2*ones(size(x,1),1);
