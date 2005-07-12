function gX = linKernDiagGradX(kern, X)

% LINKERNDIAGGRADX Gradient of linear kernel's diagonal with respect to X.

% KERN

gX = 2*X*kern.variance;
