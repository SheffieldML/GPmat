function k = linKernDiagCompute(kern, x)


% LINKERNDIAGCOMPUTE Compute diagonal of LIN kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the linear kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : linKernParamInit, kernDiagCompute, kernCreate, linKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% KERN


k =  sum(x.*x, 2)*kern.variance;
