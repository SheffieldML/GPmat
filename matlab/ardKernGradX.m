function gX = ardKernGradX(kern, x, X2)

% ARDKERNGRADX Gradient of ARD kernel with respect to a point x.

% KERN


scales = sparse(diag(kern.inputScales));

gX = kern.linearVariance.*X2*scales;

scales = sqrt(scales);

n2 = dist2(X2*scales, x*scales);
wi2 = (.5 .* kern.inverseWidth);
rbfPart = kern.rbfVariance*exp(-n2*wi2);
for i = 1:size(x, 2)
  gX(:, i) = gX(:, i) + kern.inverseWidth*kern.inputScales(i)*(X2(:, i) - x(i)).*rbfPart;
end
