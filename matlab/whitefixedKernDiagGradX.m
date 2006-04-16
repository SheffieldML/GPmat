function gX = whitefixedKernDiagGradX(kern, X)

% WHITEFIXEDKERNDIAGGRADX Gradient of white fixed noise kernel's diagonal with respect to a X.
%
% COPYRIGHT : Nathaniel J. King, 2006

% KERN

gX = whiteKernDiagGradX(kern, X);