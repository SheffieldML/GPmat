function k = translateKernDiagCompute(kern, x)

% TRANSLATEKERNDIAGCOMPUTE Compute diagonal of TRANSLATE kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the input space translation kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : translateKernParamInit, kernDiagCompute, kernCreate,
% cmpndKernDiagCompute, translateKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2007

% KERN

x = x - repmat(kern.centre, size(x, 1), 1);
k = cmpndKernDiagCompute(kern, x);
