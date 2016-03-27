function [k, sk, n2] = rbfperiodicKernCompute(kern, x, x2)

% RBFPERIODICKERNCOMPUTE Compute the RBFPERIODIC kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the RBF derived periodic
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the RBF derived periodic
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : rbfperiodicKernParamInit, kernCompute, kernCreate, rbfperiodicKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2007, 2009

% KERN
if isfield(kern, 'period')
  factor = 2*pi/kern.period;
else
  factor = 1;
end
if nargin < 3
  n2 = sin(0.5*factor*(repmat(x, 1, size(x, 1)) - repmat(x', size(x, 1), 1)));
  n2 = n2.*n2;
  wi2 = (2 .* kern.inverseWidth);
  sk = exp(-n2*wi2);
else
  n2 = sin(0.5*factor*(repmat(x, 1, size(x2, 1)) - repmat(x2', size(x, 1), 1)));  
  n2 = n2.*n2;
  wi2 = (2 .* kern.inverseWidth);
  sk = exp(-n2*wi2);
end
k = kern.variance*sk;
  
