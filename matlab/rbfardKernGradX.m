function gX = rbfardKernGradX(kern, x, X2)

% RBFARDKERNGRADX Gradient of radial basis function ARD kernel with respect to a point x.

% IVM
scales = sparse(sqrt(diag(kern.inputScales)));
gX = zeros(size(X2));
n2 = dist2(X2*scales, x*scales);
wi2 = (.5 .* kern.inverseWidth);
rbfPart = kern.variance*exp(-n2*wi2);
for i = 1:size(x, 2)
  gX(:, i) = kern.inverseWidth*kern.inputScales(i)*(X2(:, i) - x(i)).*rbfPart;
end
