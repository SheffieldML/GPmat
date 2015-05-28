function k = polyKernCompute(kern, x, x2)


% POLYKERNCOMPUTE Compute the POLY kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the polynomial
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the polynomial
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : polyKernParamInit, kernCompute, kernCreate, polyKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% KERN


if nargin < 3
  k = kern.variance*(x*x'*kern.weightVariance+kern.biasVariance).^kern.degree;
else
  k = kern.variance*(kern.weightVariance*x*x2'+kern.biasVariance).^kern.degree;
end
if issparse(x)
  k = full(k);
end
