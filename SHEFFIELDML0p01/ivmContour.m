function h = ivmContour(X, Y, Z, lineWidth)

% IVMCONTOUR Special contour plot showing decision boundary.
%
%	Description:
%
%	H = IVMCONTOUR(X, Y, Z, LINEWIDTH) plots a contour plot with a red
%	solid line at 0.5 and blue dashed lines at 0.25 and 0.75.
%	 Returns:
%	  H - handle to the contour lines.
%	 Arguments:
%	  X - input X locations (from e.g. meshgrid).
%	  Y - input Y locations.
%	  Z - input Z locations.
%	  LINEWIDTH - width of the lines to use.
%	
%	COYRIGHT : Neil D. Lawrence, 2005
%
%	See also
%	CONTOUR, IVMMESHVALS


[void, clines1] = contour(X, Y, Z, [0.25 .75], 'b--');
[void, clines2] = contour(X, Y, Z, [0.5 0.5], 'r-');
set(clines1, 'linewidth', lineWidth)
h = [clines1; clines2];
if ~isempty(clines2)
  set(clines2, 'linewidth', lineWidth)
end
