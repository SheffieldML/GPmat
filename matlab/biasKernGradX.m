function gX = biasKernGradX(kern, X, X2)

% BIASKERNGRADX Gradient of bias kernel with respect to a point x.

% KERN

gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
