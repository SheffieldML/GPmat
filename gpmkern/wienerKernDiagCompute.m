function k = wienerKernDiagCompute(kern, x)

% WIENERKERNDIAGCOMPUTE Compute diagonal of WIENER kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the wiener kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : wienerKernParamInit, kernDiagCompute, kernCreate, wienerKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2009

% KERN

k = kern.variance*x;
