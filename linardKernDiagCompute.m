function k = linardKernDiagCompute(kern, x)


% LINARDKERNDIAGCOMPUTE Compute diagonal of LINARD kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the automatic relevance determination linear kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : linardKernParamInit, kernDiagCompute, kernCreate, linardKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% KERN


scales = sparse(diag(sqrt(kern.inputScales)));
x = x*scales;

k = sum(x.*x, 2)*kern.variance;
