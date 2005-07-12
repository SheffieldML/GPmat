function gX = ardKernDiagGradX(kern, X)

% ARDKERNDIAGGRADX Gradient of ARD kernel's diagonal with respect to X.

% KERN

gX = 2*kern.linearVariance*X.*repmat(kern.inputScales, [size(X, 1) 1]);

