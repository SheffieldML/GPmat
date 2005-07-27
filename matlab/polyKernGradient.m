function g = polyKernGradient(kern, x, covGrad)

% POLYKERNGRADIENT Gradient of polynomial kernel's parameters.

% KERN

innerProd = x*x';
arg = kern.weightVariance*innerProd+kern.biasVariance;
base = kern.variance*kern.degree*arg.^(kern.degree-1);
baseCovGrad = base.*covGrad;


g(1) = sum(sum(innerProd.*baseCovGrad));
g(2) = sum(sum(baseCovGrad));
g(3) = sum(sum(covGrad.*arg.^kern.degree));

%/~
if any(isnan(g))
  warning('g is NaN')
end
%~/
