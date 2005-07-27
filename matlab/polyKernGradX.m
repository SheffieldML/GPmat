function gX = polyKernGradX(kern, X, X2)

% POLYKERNGRADX Gradient of polynomial kernel with respect to X.

% KERN

gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
for i = 1:size(X, 1);
  gX(:, :, i) = polyKernGradXpoint(kern, X(i, :), X2);
end
  

function gX = polyKernGradXpoint(kern, x, X2)

% POLYKERNGRADXPOINT Gradient with respect to one point of x.

innerProd = X2*x';  
arg = innerProd*kern.weightVariance + kern.biasVariance;
gX = zeros(size(X2));
for j = 1:size(X2, 2)
  gX(:, j) = kern.degree*kern.weightVariance*kern.variance ...
      *X2(:, j).*arg.^(kern.degree-1);
end
