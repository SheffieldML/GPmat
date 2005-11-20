function h = probitNoise3dPlot(noise, plotType, CX, CY, CZ, CZVar, varargin)

% PROBITNOISE3DPLOT Draw a 3D or contour plot for the probit.

% NOISE

fhandle = str2func(plotType);
CZ = cumGaussian((CZ + noise.bias)./sqrt(CZVar));
if nargout > 0
  h = fhandle(CX, CY, CZ, varargin{:});
else
  fhandle(CX, CY, CZ, varargin{:});
end
