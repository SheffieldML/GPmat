function probit3dPlot(noise, plotType, CX, CY, CZ, CZVar, varargin)

% PROBIT3DPLOT Draw a 3D or contour plot for the probit.

CZ = CZ + model.noise.bias;
feval(plotType, CX, CY, CZ, varargin{:});

