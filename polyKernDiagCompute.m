function k = polyKernDiagCompute(kern, x)


% POLYKERNDIAGCOMPUTE Compute diagonal of POLY kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the polynomial kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : polyKernParamInit, kernDiagCompute, kernCreate, polyKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% KERN


k =  kern.variance*(sum(x.*x, 2)*kern.weightVariance + kern.biasVariance).^kern.degree;
