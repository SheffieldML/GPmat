function h = ivmContour(X, Y, Z, lineWidth)

% IVMCONTOUR Special contour plot showing decision boundary.
% FORMAT
% DESC plots a contour plot with a red solid line at 0.5 and blue
% dashed lines at 0.25 and 0.75.
% ARG X : input X locations (from e.g. meshgrid).
% ARG Y : input Y locations.
% ARG Z : input Z locations.
% ARG lineWidth : width of the lines to use.
% RETURN h : handle to the contour lines.
% 
% SEEALSO : contour, ivmMeshVals
%
% COYRIGHT : Neil D. Lawrence, 2005

% IVM

[void, clines1] = contour(X, Y, Z, [0.25 .75], 'b--');
[void, clines2] = contour(X, Y, Z, [0.5 0.5], 'r-');
set(clines1, 'linewidth', lineWidth)
h = [clines1; clines2];
if ~isempty(clines2)
  set(clines2, 'linewidth', lineWidth)
end
