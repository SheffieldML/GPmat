function k = invcmpndKernDiagCompute(kern, x)

% INVCMPNDKERNDIAGCOMPUTE Compute diagonal of INVCMPND kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the inv. precision
% compound kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : invcmpndKernParamInit, cmpndKernDiagCompute, kernCreate, invcmpndKernCompute
%
% COPYRIGHT : Andreas C. Damianou, 2012

% KERN


% Unlike the cmpnd  version, here the nonlinearity of the inverses means
% that we cannot trivially compute the diagonal. The naive way here is to
% just compute the whole kernel and return the diagonal.

k = diag(invcmpndKernCompute(kern, x));
