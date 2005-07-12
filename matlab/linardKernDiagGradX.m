function gX = linardKernDiagGradX(kern, X)

% LINARDKERNDIAGGRADX Gradient of linear ARD kernel's diagonal with respect to X.

% KERN

gX = 2*kern.variance*X.*repmat(kern.inputScales, [size(X, 1), 1]);
