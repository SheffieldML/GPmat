function g = linKernGradient(kern, x, varargin)

% LINKERNGRADIENT Gradient of lin kernel's parameters.

% KERN

if nargin < 4
  linPart = linKernCompute(kern, x);
else
  linPart = linKernCompute(kern, x, varargin{1});
end
g(1) = sum(sum(varargin{end}.*linPart))/kern.variance;
