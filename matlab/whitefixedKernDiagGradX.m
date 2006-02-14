function gX = whitefixedKernDiagGradX(kern, X)

% WHITEFIXEDKERNDIAGGRADX Gradient of white fixed noise kernel's diagonal with respect to a X.

% KERN

gX = whiteKernDiagGradX(kern, X);