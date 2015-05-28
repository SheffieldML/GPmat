function [K, sk, n2, w2, l, l2] = gibbsKernCompute(kern, x, x2)

% GIBBSKERNCOMPUTE Compute the GIBBS kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for Mark Gibbs's non-stationary
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the Mark Gibbs's non-stationary
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : gibbsKernParamInit, kernCompute, kernCreate, gibbsKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2009

% KERN


fhandle = str2func([kern.lengthScaleTransform, 'Transform']);
l = fhandle(modelOut(kern.lengthScaleFunc, x), 'atox');
if nargin < 3
  n2 = dist2(x, x);
  L = repmat(l.*l, 1, size(l, 1));
  w2 = L + L';
  sk = ((2*l*l')./w2).^(kern.inputDimension/2).*exp(-n2./w2);
  K = kern.variance*sk;
else
  n2 = dist2(x, x2);
  l2 = fhandle(modelOut(kern.lengthScaleFunc, x2), 'atox');
  w2 = repmat(l.*l, 1, size(l2, 1))+repmat(l2.*l2, 1, size(l, 1))';
  sk = ((2*l*l2')./w2).^(kern.inputDimension/2).*exp(-n2./w2);
  K = kern.variance*sk;
end
