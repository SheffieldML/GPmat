function k = fileKernDiagCompute(kern, x)


% FILEKERNDIAGCOMPUTE Compute diagonal of FILE kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the stored file kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG index : indices of the diagonal to return.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : fileKernParamInit, kernDiagCompute, kernCreate, fileKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% KERN


k = kern.variance*fileKernRead(kern, x, 'diag');
