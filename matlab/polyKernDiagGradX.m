function gX = polyKernDiagGradX(kern, X)

% POLYKERNDIAGGRADX Gradient of polynomial kernel's diagonal with respect to X.

% KERN

gX = zeros(size(X));
for i = 1:size(X, 1);
  gX(i, :) = polyKernDiagGradXpoint(kern, X(i, :));
end
  

function gX = polyKernDiagGradXpoint(kern, x)

% POLYKERNDIAGGRADXPOINT Diagonal gradient with respect to one point of x.

innerProd = x*x';  
arg = innerProd*kern.weightVariance + kern.biasVariance;
gX = zeros(size(x));
for j = 1:size(x, 2)
  gX(:, j) = 2*kern.degree*x(:, j)*kern.weightVariance*kern.variance* ...
      arg.^(kern.degree - 1);
end
