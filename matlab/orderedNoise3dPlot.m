function h = orderedNoise3dPlot(noise, plotType, CX, CY, CZ, CZVar, varargin)

% ORDEREDNOISE3DPLOT Draw a 3D or contour plot for the probit.

% NOISE

CZ = (CZ+noise.bias)./sqrt(CZVar);
feval(plotType, CX, CY, CZ, varargin{:});
hold on
h = [];
for i = 2:noise.C-1
  CZ = CZ - noise.widths(i-1)./sqrt(CZVar);
  h = [h; feval(plotType, CX, CY, CZ, varargin{:})];
end
