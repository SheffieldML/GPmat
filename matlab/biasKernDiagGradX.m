function gX = biasKernDiagGradX(kern, X)

% BIASKERNDIAGGRADX Gradient of bias kernel's diagonal with respect to a point X.

% KERN

gX = zeros(size(X));