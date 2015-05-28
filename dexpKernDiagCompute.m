function K = dexpKernDiagCompute(kern, x)

% DEXPKERNDIAGCOMPUTE Compute diagonal of the double exponential kernel.
%
% FORMAT
% DESC computes the diagonal of the kernel matrix for the double
% exponential kernel given a column vector of inputs.
% ARG kern : the kernel structure for which the kernel matrix is computed.
% ARG x : input data in the form of a design matrix.
% RETURN K : a vector of the same size as x containing the diagonal of the
% kernel matrix computed at the given points.
%
% SEEALSO : dexpKernParamInit, kernDiagCompute, kernCreate, dexpKernCompute
%
% COPYRIGHT : David Luengo, 2009

% KERN


K = 0.5 * kern.variance * kern.decay * ones(size(x));
