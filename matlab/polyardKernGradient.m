function g = polyardKernGradient(kern, x, covGrad)

% POLYARDKERNGRADIENT Gradient of polynomial ARD kernel's parameters.

% KERN

scales = sparse(diag(sqrt(kern.inputScales)));
xScale = x*scales;
innerProd = xScale*xScale';
arg = kern.weightVariance*innerProd+kern.biasVariance;
base = kern.variance*kern.degree*arg.^(kern.degree-1);
baseCovGrad = base.*covGrad;


g(1) = sum(sum(innerProd.*baseCovGrad));
g(2) = sum(sum(baseCovGrad));
g(3) = sum(sum(covGrad.*arg.^kern.degree));

for j = 1:kern.inputDimension
  x2 = x(:, j).*x(:, j);
  g(3+j) = sum(sum((x(:, j)*x(:, j)').*baseCovGrad))*kern.weightVariance;
end


