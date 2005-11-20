function g = polyardKernGradient(kern, x, varargin)

% POLYARDKERNGRADIENT Gradient of polynomial ARD kernel's parameters.

% KERN

scales = sparse(diag(sqrt(kern.inputScales)));
xScale = x*scales;
if nargin < 4
  innerProd = xScale*xScale';
else
  xScale2 = varargin{1}*scales;
  innerProd = xScale*xScale2';
end
arg = kern.weightVariance*innerProd+kern.biasVariance;
base = kern.variance*kern.degree*arg.^(kern.degree-1);
baseCovGrad = base.*varargin{end};


g(1) = sum(sum(innerProd.*baseCovGrad));
g(2) = sum(sum(baseCovGrad));
g(3) = sum(sum(varargin{end}.*arg.^kern.degree));

if nargin < 4
  for j = 1:kern.inputDimension
    g(3+j) = sum(sum((x(:, j)*x(:, j)').*baseCovGrad))*kern.weightVariance;
  end
else
  for j = 1:kern.inputDimension
    g(3+j) = sum(sum((x(:, j)*varargin{1}(:, j)').*baseCovGrad))*kern.weightVariance;
  end
end

