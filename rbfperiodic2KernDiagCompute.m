function k = rbfperiodic2KernDiagCompute(kern, x)

% RBFPERIODIC2KERNDIAGCOMPUTE Compute diagonal of RBFPERIODIC2 kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the RBF periodic covariance with variying period kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : rbfperiodic2KernParamInit, kernDiagCompute, kernCreate, rbfperiodic2KernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2007, 2009
%
% MODIFICATIONS : Andreas C. Damianou, 2011
%
% MODIFICATIONS : Michalis K. Titsias, 2011

% KERN



k = repmat(kern.variance, size(x, 1), 1);
