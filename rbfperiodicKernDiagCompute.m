function k = rbfperiodicKernDiagCompute(kern, x)

% RBFPERIODICKERNDIAGCOMPUTE Compute diagonal of RBFPERIODIC kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the RBF derived periodic kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : rbfperiodicKernParamInit, kernDiagCompute, kernCreate, rbfperiodicKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2007

% KERN


k = repmat(kern.variance, size(x, 1), 1);
