function gX = mlpKernGradX(kern, x, x2)

% MLPKERNGRADX Gradient of multi-layer perceptron kernel with respect to a point X.

% KERN

% KERN


innerProd = x2*x';  
numer = innerProd*kern.weightVariance + kern.biasVariance;
vec1 = sum(x.*x, 2)*kern.weightVariance + kern.biasVariance + 1;
vec2 = sum(x2.*x2, 2)*kern.weightVariance + kern.biasVariance + 1;
denom = sqrt(vec2*vec1');
arg = numer./denom;
gX = zeros(size(x2));
for j = 1:size(x2, 2)
  gX(:, j)=x2(:, j)./denom...
           - vec2.*x(:, j).*numer./denom.^3;
  gX(:, j) = kern.weightVariance*kern.variance*gX(:, j)./sqrt(1-arg.*arg);
end
