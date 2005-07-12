function gX = whiteKernGradX(kern, X, X2)

% WHITEKERNGRADX Gradient of white noise kernel with respect to a point x.

% KERN

gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
