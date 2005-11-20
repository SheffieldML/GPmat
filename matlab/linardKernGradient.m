function g = linardKernGradient(kern, x, varargin)

% LINARDKERNGRADIENT Gradient of linear ARD kernel's parameters.

% KERN

g = zeros(1, size(x, 2)+1);
if nargin < 4
  k = linardKernCompute(kern, x);
else
  k = linardKernCompute(kern, x, varargin{1});
end
g(1) = sum(sum(varargin{end}.*k))/kern.variance;
if nargin < 4
  for i = 1:size(x, 2)
    g(1+i) =  x(:, i)'*varargin{end}*x(:, i)*kern.variance;
  end
else
  for i = 1:size(x, 2)
    g(1+i) =  x(:, i)'*varargin{end}*varargin{1}(:, i)*kern.variance;
  end
end