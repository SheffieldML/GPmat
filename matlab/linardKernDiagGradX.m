function gX = linardKernDiagGradX(kern, x)

% LINARDKERNDIAGGRADX Gradient of linear ARD kernel's diagonal with respect to a point x.

% IVM


gX = 2*kern.variance*x.*kern.inputScales;
