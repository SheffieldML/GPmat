function gX = ardKernDiagGradX(kern, x)

% ARDKERNDIAGGRADX Gradient of ARD kernel's diagonal with respect to a point x.

% KERN



gX = 2*kern.linearVariance*x.*kern.inputScales;

