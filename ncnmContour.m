function h = ncnmContour(X, Y, Z, lineWidth)

% NCNMCONTOUR Special contour plot showing null category region.
% FORMAT
% DESC produces a special contour plot for the null category noise
% model which places lines at the center and edges of the null
% category region.
% ARG X : input X locations for showing contours.
% ARG Y : input Y locations for showing contours.
% ARG Z : output of the null category noise model.
% ARG lineWidth : width of contour lines.
% RETURN H : handle to contour lines.

% IVM

% It's probably learnt something.
[void, clines1] =contour(X, Y, Z, [-0.5 0.5], 'b--');
[void, clines2] =contour(X, Y, Z, [0 0], 'r-');
h = [clines1; clines2];
set(h, 'linewidth', lineWidth)
