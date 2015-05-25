function k = irbfKernDiagCompute(kern, x)

% IRBFKERNDIAGCOMPUTE Compute diagonal of IRBF kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the integral of the RBF kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : irbfKernParamInit, kernDiagCompute, kernCreate, irbfKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2009

% KERN

