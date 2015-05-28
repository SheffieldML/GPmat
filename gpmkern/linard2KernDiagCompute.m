function k = linard2KernDiagCompute(kern, x)


% LINARD2KERNDIAGCOMPUTE Compute diagonal of LINARD2 kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the automatic relevance determination linear kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : linard2KernParamInit, kernDiagCompute, kernCreate, linard2KernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
%
% COPYRIGHT : Michalis K. Titsias, 2009

% KERN


scales = sparse(diag(sqrt(kern.inputScales)));
x = x*scales;

k = sum(x.*x, 2);
