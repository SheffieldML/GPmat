function [K, sk] = wienerKernCompute(kern, x, x2)

% WIENERKERNCOMPUTE Compute the WIENER kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the wiener
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the wiener
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : wienerKernParamInit, kernCompute, kernCreate, wienerKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2009

% KERN

  if any(x<0)
    error('WIENER kernel only valid for time greater than zero')
  end
  if nargin < 3
    K = repmat(x, 1, size(x, 1));
    sk = min(K, K');
  else
    if any(x2<0)
      error('WIENER kernel only valid for time greater than zero')
    end
    sk = min(repmat(x, 1, size(x2, 1)), repmat(x2', size(x, 1), 1));
  end
  K = kern.variance*sk;
end
