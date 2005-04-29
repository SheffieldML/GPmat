function h = probitNoise3dPlot(noise, plotType, CX, CY, CZ, CZVar, varargin)

% PROBITNOISE3DPLOT Draw a 3D or contour plot for the probit.

% NOISE

CZ = cumGaussian((CZ + noise.bias)./sqrt(CZVar));
if nargout > 0
  h = feval(plotType, CX, CY, CZ, varargin{:});
else
  feval(plotType, CX, CY, CZ, varargin{:});
end
