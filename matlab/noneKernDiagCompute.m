function k = noneKernDiagCompute(kern, x)


% NONEKERNDIAGCOMPUTE Compute diagonal of NONE kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the dummy kernel function kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : noneKernParamInit, kernDiagCompute, kernCreate, noneKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2008
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009

% KERN

  
k = zeros(size(x,1),1);
