function K = irbfKernCompute(kern, x, x2)

% IRBFKERNCOMPUTE Compute the IRBF kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the integral of the RBF
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the integral of the RBF
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : irbfKernParamInit, kernCompute, kernCreate, irbfKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2009

% KERN

