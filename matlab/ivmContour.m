function ivmContour(X, Y, Z, lineWidth)

% IVMCONTOUR Special contour plot showing decision boundary.

% It's probably learnt something.
[void, clines] =contour(X, Y, Z, [0.25 .75], 'b--');
set(clines, 'linewidth', lineWidth)
[void, clines] =contour(X, Y, Z, [0.5 0.5], 'r-');
set(clines, 'linewidth', lineWidth);
