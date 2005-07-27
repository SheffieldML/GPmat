function gX = polyardKernGradX(kern, X, X2)

% POLYARDKERNGRADX Gradient of polynomial ARD kernel with respect to X.

% KERN

gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
for i = 1:size(X, 1);
  gX(:, :, i) = polyardKernGradXpoint(kern, X(i, :), X2);
end
  
function gX = polyardKernGradXpoint(kern, x, X2)

% POLYARDKERNGRADXPOINT Gradient with respect to one point of x.

scales = sparse(diag(kern.inputScales));
xScaled = x*scales;
X2Scaled = X2*scales;
innerProd = X2Scaled*x';
arg = innerProd*kern.weightVariance + kern.biasVariance;
gX = zeros(size(X2));
for j = 1:size(X2, 2)
  gX(:, j) = kern.degree*kern.weightVariance*kern.inputScales(j)*kern.variance*X2(:, j).*arg.^(kern.degree-1);
end
