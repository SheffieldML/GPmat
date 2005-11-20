function g = rbfKernGradient(kern, x, varargin)

% RBFKERNGRADIENT Gradient of rbf kernel's parameters.

% KERN


% The last argument is covGrad
if nargin < 4
  [k, dist2xx] = rbfKernCompute(kern, x);
else
  [k, dist2xx] = rbfKernCompute(kern, x, varargin{1});
end
g(1) = - .5*sum(sum(varargin{end}.*k.*dist2xx));
g(2) =  sum(sum(varargin{end}.*k))/kern.variance;
%/~
if any(isnan(g))
  warning('g is NaN')
end
%~/