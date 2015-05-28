function [k, sk, n2] = rbfperiodic2KernCompute(kern, x, x2)

% RBFPERIODIC2KERNCOMPUTE Compute the RBFPERIODIC2 kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the RBF periodic covariance with variying period
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the RBF periodic covariance with variying period
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : rbfperiodic2KernParamInit, kernCompute, kernCreate, rbfperiodic2KernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2007, 2009
%
% MODIFICATIONS : Andreas C. Damianou, 2011
%
% MODIFICATIONS : Michalis K. Titsias, 2011

% KERN

factor = kern.factor; % Default (if period is fixed: 2*pi/kern.period)
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
% Test kernel with: kernTest('rbfperiodic2',1)

  
