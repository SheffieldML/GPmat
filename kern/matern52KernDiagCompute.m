function k = matern52KernDiagCompute(kern, x)

% MATERN52KERNDIAGCOMPUTE Compute diagonal of MATERN52 kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the matern kernel with nu=5/2 kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : matern52KernParamInit, kernDiagCompute, kernCreate, matern52KernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

k = repmat(kern.variance, size(x, 1), 1);
