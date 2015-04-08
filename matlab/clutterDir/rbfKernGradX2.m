function gX = rbfKernGradX2(kern, X, X2)

% RBFKERNGRADX Gradient of Radial basis function kernel with respect to X.


gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
for i = 1:size(X, 1);
  gX(:, :, i) = rbfKernGradXpoint(kern, X(i, :), X2);
end
  

function gX = rbfKernGradXpoint(kern, x, X2)

% RBFKERNGRADXPOINT Gradient with respect to one point of x.

gX = zeros(size(X2));
n2 = dist2(X2, x);
wi2 = (.5 .* kern.inverseWidth);
rbfPart = kern.variance*exp(-n2*wi2);
for i = 1:size(x, 2)
  gX(:, i) = kern.inverseWidth*(X2(:, i) - x(i)).*rbfPart;
end
