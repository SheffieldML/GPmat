function gX = biasKernDiagGradX(kern, x)

% BIASKERNDIAGGRADX Gradient of bias kernel's diagonal with respect to a point x.

% IVM

gX = zeros(size(x));