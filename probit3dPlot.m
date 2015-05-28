function probit3dPlot(noise, plotType, CX, CY, CZ, CZVar, varargin)

% PROBIT3DPLOT Draw a 3D or contour plot for the probit.

% NOISE

CZ = CZ + model.noise.bias;
fhandle = str2func(plotType);
fhandle(CX, CY, CZ, varargin{:});

