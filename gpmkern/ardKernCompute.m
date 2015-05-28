function [k, rbfPart, linearPart, n2] = ardKernCompute(kern, x, x2)


% ARDKERNCOMPUTE Compute the ARD kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the pre-built RBF and linear ARD
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the pre-built RBF and linear ARD
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : ardKernParamInit, kernCompute, kernCreate, ardKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2004

% KERN


scales = diag(sqrt(kern.inputScales));
x = x*scales;
    
if nargin < 3
  n2 = dist2(x, x);
  wi2 = (.5 .* kern.inverseWidth);
  rbfPart = kern.rbfVariance*exp(-n2*wi2);
  linearPart = x*x'*kern.linearVariance;
  k = rbfPart + kern.whiteVariance*eye(size(x, 1)) + linearPart;
else
  x2 = x2*scales;
  n2 = dist2(x, x2);
  wi2 = (.5 .* kern.inverseWidth);
  rbfPart = kern.rbfVariance*exp(-n2*wi2);
  linearPart = x*x2'*kern.linearVariance;
  k = rbfPart + linearPart;
end
k = k + kern.biasVariance;
