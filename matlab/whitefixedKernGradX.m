function gX = whitefixedKernGradX(kern, X, X2)

% WHITEFIXEDKERNGRADX Gradient of white fixed noise kernel with respect to a point x.
%
% COPYRIGHT : Nathaniel J. King, 2006

% KERN

gX = whiteKernGradX(kern, X, X2);