function gX = rbfKernGradX(kern, x, x2)

% RBFKERNGRADX Gradient of Radial basis function kernel with respect to a point X.

% IVM

gX = zeros(size(x2));
n2 = dist2(x2, x);
wi2 = (.5 .* kern.inverseWidth);
rbfPart = kern.variance*exp(-n2*wi2);
for i = 1:size(x, 2)
  gX(:, i) = kern.inverseWidth*(x2(:, i) - x(i)).*rbfPart;
end
