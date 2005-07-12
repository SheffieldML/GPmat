function gX = mlpardKernDiagGradX(kern, X)

% MLPARDKERNDIAGGRADX Gradient of multi-layer perceptron ARD kernel's diagonal with respect to X.

% KERN

gX = zeros(size(X));
for i = 1:size(X, 1);
  gX(i, :) = mlpardKernDiagGradXpoint(kern, X(i, :));
end
  

function gX = mlpardKernDiagGradXpoint(kern, x)

% MLPARDKERNDIAGGRADXPOINT Diagonal gradient with respect to one point of x.

innerProd = x*sparse(diag(kern.inputScales))*x';  
numer = innerProd*kern.weightVariance + kern.biasVariance;
denom = numer + 1;
arg = numer./denom;
gX = zeros(size(x));
for j = 1:size(x, 2)
  gX(:, j)=1./denom...
           - numer./denom.^2;
  gX(:, j) = 2*kern.inputScales(j)*x(:, j)*kern.weightVariance*kern.variance*gX(:, j)./sqrt(1-arg.*arg);
end
