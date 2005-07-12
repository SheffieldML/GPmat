function gX = ardKernGradX(kern, X, X2)

% ARDKERNGRADX Gradient of ARD kernel with respect to X.

% KERN

gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
for i = 1:size(X, 1);
  gX(:, :, i) = ardKernGradXpoint(kern, X(i, :), X2);
end

function gX = ardKernGradXpoint(kern, x, X2)

% ARDKERNGRADXPOINT Gradient with respect to one point of x.

scales = sparse(diag(kern.inputScales));

gX = kern.linearVariance.*X2*scales;

scales = sqrt(scales);

n2 = dist2(X2*scales, x*scales);
wi2 = (.5 .* kern.inverseWidth);
rbfPart = kern.rbfVariance*exp(-n2*wi2);
for i = 1:size(x, 2)
  gX(:, i) = gX(:, i) + kern.inverseWidth*kern.inputScales(i)*(X2(:, i) - x(i)).*rbfPart;
end
