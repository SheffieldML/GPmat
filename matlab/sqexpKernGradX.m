function gX = sqexpKernGradX(kern, x, x2)

% SQEXPKERNGRADX Gradient of squared exponential kernel with respect to a point X.

% IVM

gX = zeros(size(x2));
n2 = dist2(x2, x);
wi2 = (.5 .* kern.inverseWidth);
rbfPart = kern.rbfVariance*exp(-n2*wi2);
for i = 1:size(x, 2)
  gX(:, i) = kern.inverseWidth*(x2(:, i) - x(i)).*rbfPart;
end
