function gX = mlpKernGradX(kern, X, X2)

% MLPKERNGRADX Gradient of multi-layer perceptron kernel with respect to X.

% KERN

gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
for i = 1:size(X, 1);
  gX(:, :, i) = mlpKernGradXpoint(kern, X(i, :), X2);
end
  

function gX = mlpKernGradXpoint(kern, x, X2)

% MLPKERNGRADXPOINT Gradient with respect to one point of x.

innerProd = X2*x';  
numer = innerProd*kern.weightVariance + kern.biasVariance;
vec1 = sum(x.*x, 2)*kern.weightVariance + kern.biasVariance + 1;
vec2 = sum(X2.*X2, 2)*kern.weightVariance + kern.biasVariance + 1;
denom = sqrt(vec2*vec1');
arg = numer./denom;
gX = zeros(size(X2));
for j = 1:size(X2, 2)
  gX(:, j)=X2(:, j)./denom...
           - vec2.*x(:, j).*numer./denom.^3;
  gX(:, j) = kern.weightVariance*kern.variance*gX(:, j)./sqrt(1-arg.*arg);
end
