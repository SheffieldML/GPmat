function [h1, h2] = ivm3dPlot(model, plotType, iter, X, y)

% IVM3DPLOT Make a 3-D or contour plot of the IVM.

% IVM

h2 = [];
lineWidth = 2;
if nargin < 5
  X = model.X;
  y = model.y;
end
clf
h1 = noisePointPlot(model.noise, X, y, 'times', 20, 10, lineWidth);
title(['Iteration ' num2str(iter)])

if length(model.I)>0
  % It's probably learnt something.
  xlim = get(gca, 'xlim');
  ylim = get(gca, 'ylim');
  [CX, CY, CZ, CZVar] = ivmMeshVals(model, xlim, ylim, 60);
  switch plotType
   case {'ivmContour', 'ncnmContour'}
    threeDargs = cell(1);
    threeDargs{1} = 2;
   otherwise
    threeDargs = cell(0);
  end
  h2 =noise3dPlot(model.noise, plotType, CX, CY, CZ, CZVar, threeDargs{:});   
end
drawnow