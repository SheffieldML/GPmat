function gX = whitefixedKernGradX(kern, X, X2)

% WHITEFIXEDKERNGRADX Gradient of white fixed noise kernel with respect to a point x.

% KERN

gX = whiteKernGradX(kern, X, X2);