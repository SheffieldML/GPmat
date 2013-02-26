function h = ncnmNoise3dPlot(noise, plotType, CX, CY, CZ, CZVar, varargin)

% NCNMNOISE3DPLOT Draw a 3D or contour plot for the NCNM noise model.
%
%	Description:
%	h = ncnmNoise3dPlot(noise, plotType, CX, CY, CZ, CZVar, varargin)
%

if nargout > 0
  h = feval(plotType, CX, CY, CZ, varargin{:});
else
  feval(plotType, CX, CY, CZ, varargin{:});
end

