function gX = rbfardKernDiagGradX(kern, X)

% RBFARDKERNDIAGGRADX Gradient of radial basis function ARD kernel's diagonal with respect to X.

% KERN

gX = zeros(size(X));