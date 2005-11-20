function g = biasKernGradient(kern, x, varargin)

% BIASKERNGRADIENT Gradient of bias kernel's parameters.

% KERN

g = sum(sum(varargin{end}));
