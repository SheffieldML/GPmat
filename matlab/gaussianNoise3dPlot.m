function gaussianNoise3dPlot(noise, plotType, CX, CY, CZ, CZVar, varargin)

% GAUSSIANNOISE3DPLOT Draw a 3D or contour plot for the Gassian noise model.

% NOISE

CZ = (CZ+noise.bias);
feval(plotType, CX, CY, CZ, varargin{:});
