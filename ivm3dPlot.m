function [h1, h2] = ivm3dPlot(model, plotType, iter, X, y)

% IVM3DPLOT Make a 3-D or contour plot of the IVM.
% FORMAT
% DESC makes a 3-D or contour plot from an IVM model to show
% results.
% ARG model : the model from which the plot is generated.
% ARG plotType : type of plot to be used, options include
% 'ncnmContour', which gives the contours at the edge of the null
% category region for the null category noise model, and
% 'ivmContour'.
% ARG iter : iteration number, if give it is used at the title of
% the plot.
% ARG X : optional argument, if given it is the plotted input
% locations, otherwise model.X is used.
% ARG y : optional argument, if given it is the plotted target
% locations, otherwise model.y is used.
% RETURN h1 : a handle to the data in the plot.
% RETURN h2 : a handle to any contours, or functions plotted
% derived from the IVM.
%
% SEEALSO : noise3dPlot, noisePointPlot, ivmMeshVals, ivmCreate
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% IVM

h2 = [];
lineWidth = 2;
if nargin < 5
  X = model.X;
  y = model.y;
end
if isempty(X)
  X = model.X;
end
if isempty(y)
  y = model.y;
end
clf
h1 = noisePointPlot(model.noise, X, y, 'times', 20, 10, lineWidth);

if ~isempty(iter)
  title(['Iteration ' num2str(iter)])
end

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
  h2 = noise3dPlot(model.noise, plotType, CX, CY, CZ, CZVar, threeDargs{:});   
end
drawnow
