function probitNoise3dPlot(noise, plotType, CX, CY, CZ, CZVar, varargin)

% PROBITNOISE3DPLOT Draw a 3D or contour plot for the probit.

CZ = cumGaussian((CZ + noise.bias)./sqrt(CZVar));
feval(plotType, CX, CY, CZ, varargin{:});

