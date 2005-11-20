function g = rbfardKernGradient(kern, x, varargin)

% RBFARDKERNGRADIENT Gradient of radial basis function ARD kernel's parameters.

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