function k = whitehKernDiagCompute(kern, x)

% WHITEHKERNDIAGCOMPUTE Compute diagonal of WHITEH kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the whiteh noise kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : whitehKernParamInit, kernDiagCompute, kernCreate, whitehKernCompute
%
% COPYRIGHT : Mauricio A. Alvarez, Neil D. Lawrence, 2009

% KERN

% /~MAURICIO : Intended for the School Data only
% ~/
k = kern.variance./x;

