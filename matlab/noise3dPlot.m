function noise3dPlot(noise, plotType, CX, CY, CZ, CZVar, varargin)

% NOISE3DPLOT Draw a 3D or contour plot for the relevant noise model.

functionName = [noise.type 'Noise3dPlot'];
if exist(functionName) == 2
  feval(functionName, noise, plotType, CX, CY, CZ, CZVar, varargin{:});
end

