function gX = polyardKernDiagGradX(kern, X)

% POLYARDKERNDIAGGRADX Gradient of polynomial ARD kernel's diagonal with respect to X.

% KERN

gX = zeros(size(X));
for i = 1:size(X, 1);
  gX(i, :) = polyardKernDiagGradXpoint(kern, X(i, :));
end
  

function gX = polyardKernDiagGradXpoint(kern, x)

% POLYARDKERNDIAGGRADXPOINT Diagonal gradient with respect to one point of x.

innerProd = x*sparse(diag(kern.inputScales))*x';  
arg = innerProd*kern.weightVariance + kern.biasVariance;
gX = zeros(size(x));
for j = 1:size(x, 2)
  gX(:, j) = kern.degree*2*kern.inputScales(j)*x(:, j)*kern.weightVariance*kern.variance.*arg.^(kern.degree-1);
end
