function g = fileKernGradient(kern, x, varargin)

% FILEKERNGRADIENT Gradient of file stored kernel's parameters.

% KERN

% The last argument is covGrad
if nargin < 4
  k = fileKernCompute(kern, x);
else
  k = fileKernCompute(kern, x, varargin{1});
end
g(1) =  sum(sum(varargin{end}.*k))/kern.variance;
%/~
if any(isnan(g))
  warning('g is NaN')
end
%~/