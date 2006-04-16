function k = ardKernDiagCompute(kern, x)


% ARDKERNDIAGCOMPUTE Compute diagonal of ARD kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the pre-built RBF and linear ARD kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : ardKernParamInit, kernDiagCompute, kernCreate, ardKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2004

% KERN


scales = sparse(diag(sqrt(kern.inputScales)));
x = x*scales;

rbfPart = ones(size(x, 1), 1);
linearPart = sum(x.*x, 2)*kern.linearVariance;
k = rbfPart*(kern.rbfVariance + kern.whiteVariance) + kern.biasVariance +linearPart;
