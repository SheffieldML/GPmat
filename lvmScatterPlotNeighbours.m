function [ax, data] = lvmScatterPlotNeighbours(model, YLbls, ax);

% LVMSCATTERPLOTNEIGHBOURS 2-D scatter plot of the latent points with neighbourhood.
% FORMAT
% DESC produces a visualisation of the latent space with the given model.
% ARG model : the model for which the scatter plot is being produced.
% RETURN ax : the axes handle where the scatter plot was placed.
%
% DESC produces a visualisation of the latent space for the given model, 
% using the provided labels to distinguish the latent points.
% ARG model : the model for which the scatter plot is being produced.
% ARG lbls : labels for each data point so that they may be given different
% symbols. Useful when each data point is associated with a different
% class.
% RETURN ax : the axes handle where the scatter plot was placed.
% 
% DESC produces a visualisation of the latent space for the given model, 
% using the provided labels to distinguish the latent points.
% ARG model : the model for which the scatter plot is being produced.
% ARG lbls : labels for each data point so that they may be given different
% symbols. Useful when each data point is associated with a different
% class.
% ARG ax : the axes where the plot is to be placed.
% RETURN ax : the axes handle where the scatter plot was placed.
%
% COPYRIGHT : Neil D. Lawrence, 2008
%
% SEEALSO : lvmScatterPlot, lvmTwoDPlot, lvmScatterPlotColor
 
% MLTOOLS

if nargin<3
  ax = [];
  if nargin < 2
    YLbls = [];
  end
end

if isempty(YLbls)
  symbol = [];
else
  symbol = getSymbols(size(YLbls,2));
end

[ax, data] = lvmScatterPlot(model, YLbls, ax);


for i = 1:model.N
  for j = 1:length(model.indices(i, :))
    line([model.X(i, 1)  ...
         model.X(model.indices(i, j), 1)], ...
         [model.X(i, 2) ...
         model.X(model.indices(i, j), 2)]);
  end
end
