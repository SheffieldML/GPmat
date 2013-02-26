function h = ncnmContour(X, Y, Z, lineWidth)

% NCNMCONTOUR Special contour plot showing null category region.
%
%	Description:
%	h = ncnmContour(X, Y, Z, lineWidth)
%

% It's probably learnt something.
[void, clines1] =contour(X, Y, Z, [-0.5 0.5], 'b--');
[void, clines2] =contour(X, Y, Z, [0 0], 'r-');
h = [clines1; clines2];
set(h, 'linewidth', lineWidth)
