function gX = rbfardKernGradX(kern, X, X2)

% RBFARDKERNGRADX Gradient of radial basis function ARD kernel with respect to a X.

% KERN

gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
for i = 1:size(X, 1);
  gX(:, :, i) = rbfardKernGradXpoint(kern, X(i, :), X2);
end
  

function gX = rbfardKernGradXpoint(kern, x, X2)

% RBFARDKERNGRADXPOINT Gradient with respect to one point of x.

scales = sparse(sqrt(diag(kern.inputScales)));
gX = zeros(size(X2));
n2 = dist2(X2*scales, x*scales);
wi2 = (.5 .* kern.inverseWidth);
rbfPart = kern.variance*exp(-n2*wi2);
for i = 1:size(x, 2)
  gX(:, i) = kern.inverseWidth*kern.inputScales(i)*(X2(:, i) - x(i)).*rbfPart;
end
