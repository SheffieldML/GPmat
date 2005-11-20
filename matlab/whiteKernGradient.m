function g = whiteKernGradient(kern, x, varargin)

% WHITEKERNGRADIENT Gradient of white noise kernel's parameters.

% KERN

if nargin < 4
  g = trace(varargin{end});
else
  g = 0;
end
