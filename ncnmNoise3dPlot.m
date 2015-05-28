function h = ncnmNoise3dPlot(noise, plotType, CX, CY, CZ, CZVar, varargin)

% NCNMNOISE3DPLOT Draw a 3D or contour plot for the NCNM noise model.
% FORMAT
% DESC draws a contour plot of outputs for the null category noise
% model.
% ARG noise : structure containing the noise model.
% ARG plotType : the type of plot to create (typically 'ncnmContour').
% ARG X : input X locations for showing contours.
% ARG Y : input Y locations for showing contours.
% ARG Z : input means to the null category noise model.
% ARG CZVar : input variances to the null category noise model.
% ARG P1, P2, ... : optional input arguments to the plot type.
% RETURN H : handle to contour lines.

% NOISE

if nargout > 0
  h = feval(plotType, CX, CY, CZ, varargin{:});
else
  feval(plotType, CX, CY, CZ, varargin{:});
end

