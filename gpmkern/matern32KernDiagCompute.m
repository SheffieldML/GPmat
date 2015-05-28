function k = matern32KernDiagCompute(kern, x)

% MATERN32KERNDIAGCOMPUTE Compute diagonal of MATERN32 kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the matern kernel with nu=3/2 kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : matern32KernParamInit, kernDiagCompute, kernCreate, matern32KernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

k = repmat(kern.variance, size(x, 1), 1);
