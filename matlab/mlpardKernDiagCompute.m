function k = mlpardKernDiagCompute(kern, x)


% MLPARDKERNDIAGCOMPUTE Compute diagonal of MLPARD kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the automatic relevance determination multi-layer perceptron kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : mlpardKernParamInit, kernDiagCompute, kernCreate, mlpardKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% KERN


scales = sparse(diag(sqrt(kern.inputScales)));
x = x*scales;
numer = sum(x.*x, 2)*kern.weightVariance + kern.biasVariance;
denom = numer+1;
k = kern.variance*2/pi*asin(numer./denom);
