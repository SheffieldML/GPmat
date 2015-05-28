function g = rbfardKernGradient(kern, x, varargin)

% RBFARDKERNGRADIENT Gradient of RBFARD kernel's parameters.
% FORMAT
% DESC computes the gradient of functions with respect to the
% automatic relevance determination radial basis function
% kernel's parameters. As well as the kernel structure and the
% input positions, the user provides a matrix PARTIAL which gives
% the partial derivatives of the function with respect to the
% relevant elements of the kernel matrix. 
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG x : the input locations for which the gradients are being
% computed. 
% ARG partial : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The argument takes
% the form of a square matrix of dimension  numData, where numData is
% the number of rows in X.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters. The ordering of the vector should match
% that provided by the function kernExtractParam.
%
% FORMAT
% DESC computes the derivatives as above, but input locations are
% now provided in two matrices associated with rows and columns of
% the kernel matrix. 
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG x1 : the input locations associated with the rows of the
% kernel matrix.
% ARG x2 : the input locations associated with the columns of the
% kernel matrix.
% ARG partial : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The matrix should
% have the same number of rows as X1 and the same number of columns
% as X2 has rows.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters.
%
% SEEALSO rbfardKernParamInit, kernGradient, rbfardKernDiagGradient, kernGradX
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% KERN


g = zeros(1, size(x, 2)+2);

if nargin < 4
  [k, dist2xx] = rbfardKernCompute(kern, x);
else
  [k, dist2xx] = rbfardKernCompute(kern, x, varargin{1});
end
covGradK = varargin{end}.*k;
g(1) = - .5*sum(sum(covGradK.*dist2xx));
g(2) =  sum(sum(varargin{end}.*k))/kern.variance;

if nargin < 4
  for i = 1:size(x, 2)
    g(2 + i)  =  -(sum(covGradK*(x(:, i).*x(:, i))) ...
                   -x(:, i)'*covGradK*x(:, i))*kern.inverseWidth;
  end
else
  for i = 1:size(x, 2)
    g(2 + i) = -(0.5*sum(covGradK'*(x(:, i).*x(:, i))) ...
                 + 0.5*sum(covGradK*(varargin{1}(:, i).*varargin{1}(:, i)))...
                 -x(:, i)'*covGradK*varargin{1}(:, i))*kern.inverseWidth;
  end
end
