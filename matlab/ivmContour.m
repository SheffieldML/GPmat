function h = ivmContour(X, Y, Z, lineWidth)

% IVMCONTOUR Special contour plot showing decision boundary.

% IVM

% It's probably learnt something.
[void, clines1] =contour(X, Y, Z, [0.25 .75], 'b--');
[void, clines2] =contour(X, Y, Z, [0.5 0.5], 'r-');
h = [clines1; clines2];
set(h, 'linewidth', lineWidth)
