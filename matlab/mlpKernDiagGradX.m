function gX = mlpKernDiagGradX(kern, X)

% MLPKERNDIAGGRADX Gradient of multi-layer perceptron kernel's diagonal with respect to X.

% KERN

gX = zeros(size(X));
for i = 1:size(X, 1);
  gX(i, :) = mlpKernDiagGradXpoint(kern, X(i, :));
end
  

function gX = mlpKernDiagGradXpoint(kern, x)

% MLPKERNDIAGGRADXPOINT Diagonal gradient with respect to one point of x.

innerProd = x*x';  
numer = innerProd*kern.weightVariance + kern.biasVariance;
denom = numer + 1;
arg = numer./denom;
gX = zeros(size(x));
for j = 1:size(x, 2)
  gX(:, j)=1./denom...
           - numer./denom.^2;
  gX(:, j) = 2*x(:, j)*kern.weightVariance*kern.variance*gX(:, j)./sqrt(1-arg.*arg);
end
