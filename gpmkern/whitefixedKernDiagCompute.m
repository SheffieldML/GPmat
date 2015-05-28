function k = whitefixedKernDiagCompute(kern, x)


% WHITEFIXEDKERNDIAGCOMPUTE Compute diagonal of WHITEFIXED kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the fixed parameter white noise kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : whitefixedKernParamInit, kernDiagCompute, kernCreate, whitefixedKernCompute
%
% COPYRIGHT : Nathaniel J. King, 2006

% KERN


k = whiteKernDiagCompute(kern, x);
