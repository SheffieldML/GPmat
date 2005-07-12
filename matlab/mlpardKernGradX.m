function gX = mlpardKernGradX(kern, X, X2)

% MLPARDKERNGRADX Gradient of multi-layer perceptron ARD kernel with respect to a X.

% KERN

gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
for i = 1:size(X, 1);
  gX(:, :, i) = mlpardKernGradXpoint(kern, X(i, :), X2);
end
  
function gX = mlpardKernGradXpoint(kern, x, X2)

% MLPARDKERNGRADXPOINT Gradient with respect to one point of x.

scales = sparse(diag(kern.inputScales));
xScaled = x*scales;
X2Scaled = X2*scales;
innerProd = X2Scaled*x';
numer = innerProd*kern.weightVariance + kern.biasVariance;
vec1 = sum(xScaled.*x, 2)*kern.weightVariance + kern.biasVariance + 1;
vec2 = sum(X2Scaled.*X2, 2)*kern.weightVariance + kern.biasVariance + 1;
denom = sqrt(vec2*vec1');
arg = numer./denom;
gX = zeros(size(X2));
for j = 1:size(X2, 2)
  gX(:, j)=X2(:, j)./denom...
           - vec2.*x(:, j).*numer./denom.^3;
  gX(:, j) = kern.weightVariance*kern.inputScales(j)*kern.variance*gX(:, j)./sqrt(1-arg.*arg);
end
