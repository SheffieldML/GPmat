function gX = rbfardKernDiagGradX(kern, x)

% RBFARDKERNDIAGGRADX Gradient of radial basis function ARD kernel's diagonal with respect to a point x.

% IVM

gX = zeros(size(x));