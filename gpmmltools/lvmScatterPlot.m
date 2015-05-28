function [ax, data] = lvmScatterPlot(model, YLbls, ax, dims, defaultVals);

% LVMSCATTERPLOT 2-D scatter plot of the latent points.
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
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
%
% SEEALSO : lvmVisualise, lvmTwoDPlot, lvmScatterPlotColor
 
% MLTOOLS
  if nargin < 5
    defaultVals = zeros(1, size(model.X, 2));
    
    if nargin < 4
      dims = [1, 2];
      if nargin<3
        ax = [];
        if nargin < 2
          YLbls = [];
        end
      end
    end
  end
  if isempty(YLbls)
    symbol = getSymbols(1);
  else
    if iscell(YLbls)
      symbol = getSymbols(size(YLbls{1},2));
    else
      symbol = getSymbols(size(YLbls,2));
    end
  end
  x1Min = min(model.X(:, dims(1)));
  x1Max = max(model.X(:, dims(1)));
  x1Span = x1Max - x1Min;
  x1Min = x1Min - 0.05*x1Span;
  x1Max = x1Max + 0.05*x1Span;
  x1 = linspace(x1Min, x1Max, 150);
  
  x2Min = min(model.X(:, dims(2)));
  x2Max = max(model.X(:, dims(2)));
  x2Span = x2Max - x2Min;
  x2Min = x2Min - 0.05*x2Span;
  x2Max = x2Max + 0.05*x2Span;
  x2 = linspace(x2Min, x2Max, 150);
  
  %if size(model.X, 2)==2
    
    [X1, X2] = meshgrid(x1, x2);
    XTest = repmat(defaultVals, prod(size(X1)), 1);
    XTest(:, dims(1)) = X1(:);
    XTest(:, dims(2)) = X2(:);
    varsigma = modelPosteriorVar(model, XTest);
  try
    [X1, X2] = meshgrid(x1, x2);
    XTest = repmat(defaultVals, prod(size(X1)), 1);
    XTest(:, dims(1)) = X1(:);
    XTest(:, dims(2)) = X2(:);
    varsigma = modelPosteriorVar(model, XTest);  
    posteriorVarDefined = true;
  catch 
    [lastMsg, lastId] = lasterr;
    disp(lastId)
    if isoctave || strcmp(lastId, 'MATLAB:UndefinedFunction')
      posteriorVarDefined = false;
    else
      rethrow(lasterror);
    end
  end
  if posteriorVarDefined
    d = model.d;
    if size(varsigma, 2) == 1
      dataMaxProb = -0.5*d*log(varsigma);
    else
      dataMaxProb = -.5*sum(log(varsigma), 2);
    end
    
    if isempty(ax)
      figure(1)
      clf
      % Create the plot for the data
      ax = axes('position', [0.05 0.05 0.9 0.9]);
    else
      axes(ax);
    end
    hold on
    
    C = reshape(dataMaxProb, size(X1));
    
    % Rescale it
    C = C - min(min(C));
    if max(max(C))~=0
      C = C/max(max(C));
      C = round(C*63);
      image(x1, x2, C);
    end
    
    %[c, h] = contourf(X1, X2, log10(reshape(1./varsigma(:, 1), size(X1))), 128); 
    % shading flat
    colormap gray;
    %colorbar
  end
  
  data = lvmTwoDPlot(model.X(:, dims), YLbls, symbol);
  switch model.type
   case 'dnet'
    plot(model.X_u(:, dims(1)), model.X_u(:, dims(2)), 'g.')
  end
% elseif size(model.X, 2)==3
%   x3Min = min(model.X(:, 3));
%   x3Max = max(model.X(:, 3));
%   x3Span = x3Max - x3Min;
%   x3Min = x3Min - 0.05*x3Span;
%   x3Max = x3Max + 0.05*x3Span;
%   x3 = linspace(x3Min, x3Max, 150);

%   data = lvmThreeDPlot(model.X, YLbls, symbol);
% end

xLim = [min(x1) max(x1)];
yLim = [min(x2) max(x2)];
set(gca, 'xlim', xLim);
set(gca, 'ylim', yLim);
% if size(model.X, 2) == 3
%   zLim = [min(x3) max(x3)];
%   set(ax, 'zLim', zLim);
% end
set(gca, 'fontname', 'arial');
set(gca, 'fontsize', 20);

ax = gca;
