function g = polyKernGradient(kern, x, varargin)

% POLYKERNGRADIENT Gradient of polynomial kernel's parameters.

% KERN
if nargin < 4
  innerProd = x*x';
else
  innerProd = x*varargin{1}';
end
arg = kern.weightVariance*innerProd+kern.biasVariance;
base = kern.variance*kern.degree*arg.^(kern.degree-1);
baseCovGrad = base.*varargin{end};


g(1) = sum(sum(innerProd.*baseCovGrad));
g(2) = sum(sum(baseCovGrad));
g(3) = sum(sum(varargin{end}.*arg.^kern.degree));

%/~
if any(isnan(g))
  warning('g is NaN')
end
%~/
